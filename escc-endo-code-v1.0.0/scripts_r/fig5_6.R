# Figure 5 & 6: clinical association and diagnostic panel
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(readxl); library(scales)
})
FIG <- "D:/escc_sc/figs"; DATA <- "D:/escc_sc/data"
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 9, face = "bold"))

# ---------------- Fig 5: proportion vs stage ----------------
cl <- read_excel(file.path(DATA, "GSE160269_SupplData1_clinical.xlsx"), sheet = "S1a_scCohort", skip = 1)
names(cl)[1] <- "patient"
cl$patient <- trimws(as.character(cl$patient))
cl$stage <- factor(cl[["Pathologic stage"]], levels = c("I", "II", "III"))
cl$Tstage <- as.numeric(gsub("[^0-9]", "", cl[["T stage"]]))
cl$Nstage <- as.numeric(gsub("[^0-9]", "", cl[["N stage"]]))

disp <- c(Tip_angiogenic = "Tip-like", IFNgamma_EC = "Interferon-responsive",
          Matrix_EC = "Matrix-remodeling", Venous_inflamed = "Inflamed venule",
          Venous = "Venous", Capillary = "Capillary", Arterial = "Arterial")
prop <- read.csv(file.path(DATA, "celltype_proportion_closed.csv")) %>%
  filter(!grepl("^Contam", celltype)) %>% filter(type == "Tumor") %>%
  mutate(patient = sub("T$", "", sample),
         celltype = ifelse(celltype %in% names(disp),
                           unname(disp[as.character(celltype)]), as.character(celltype)))

key <- unname(disp[c("Tip_angiogenic", "IFNgamma_EC", "Matrix_EC",
                     "Venous_inflamed", "Venous", "Capillary", "Arterial")])
# 未检出的状态记为 0
grid <- expand.grid(sample = unique(prop$sample), celltype = key, stringsAsFactors = FALSE)
prop0 <- grid %>%
  left_join(prop[, c("sample", "celltype", "n", "frac")], by = c("sample", "celltype")) %>%
  mutate(frac = ifelse(is.na(frac), 0, frac), n = ifelse(is.na(n), 0, n),
         patient = sub("T$", "", sample))
m <- prop0 %>% inner_join(cl, by = "patient")
d <- m %>% filter(celltype %in% key) %>% mutate(celltype = factor(celltype, levels = key))

lab <- d %>% group_by(celltype) %>%
  summarise(rho = suppressWarnings(cor(frac, as.numeric(stage), method = "spearman", use = "complete.obs")),
            p = suppressWarnings(cor.test(frac, as.numeric(stage), method = "spearman")$p.value), .groups = "drop") %>%
  mutate(txt = sprintf("\u03c1 = %.2f\nP = %.2g", rho, p))

lab$q_stage_BH <- p.adjust(lab$p, method = "BH")
write.csv(lab, "D:/escc_sc/out_v37/stage_lab_7states.csv", row.names = FALSE)

p5a <- ggplot(d, aes(stage, frac, fill = celltype)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.25, width = 0.66) +
  scale_fill_manual(values = c("#b03a2e","#7b4fa0","#d98b5f","#c6a0d9","#8ec6d9","#4a9bc4","#1b6ca8")) +
  scale_y_continuous(labels = percent) +
  facet_wrap(~celltype, nrow = 1, scales = "free_y",
             labeller = labeller(celltype = c(
               "Tip-like" = "Tip-like",
               "Interferon-responsive" = "IFN-\nresponsive",
               "Matrix-remodeling" = "Matrix-\nremodeling",
               "Inflamed venule" = "Inflamed\nvenule",
               "Venous" = "Venous",
               "Capillary" = "Capillary",
               "Arterial" = "Arterial"))) +
  geom_text(data = lab, aes(x = 2, y = Inf, label = txt), vjust = 1.1, size = 1.9, inherit.aes = FALSE) +
  labs(x = "Pathological stage", y = "Fraction of endothelial cells") +
  th + theme(legend.position = "none", strip.background = element_blank(),
             strip.text = element_text(size = 6.0, lineheight = 0.9),
             axis.text.x = element_text(size = 6.5))

# tumour grade / T / N summary heat
summ <- d %>% group_by(celltype) %>%
  summarise(stage = cor(frac, as.numeric(stage), method = "spearman", use = "complete.obs"),
            T = cor(frac, Tstage, method = "spearman", use = "complete.obs"),
            N = cor(frac, Nstage, method = "spearman", use = "complete.obs"), .groups = "drop") %>%
  pivot_longer(-celltype, names_to = "var", values_to = "rho")
write.csv(summ, "D:/escc_sc/out_v37/stage_summ_7states.csv", row.names = FALSE)
p5b <- ggplot(summ, aes(var, celltype, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2) +
  scale_fill_gradient2(low = "#1b6ca8", mid = "white", high = "#b03a2e", midpoint = 0,
                       limits = c(-0.5, 0.5)) +
  labs(x = NULL, y = NULL, fill = "Spearman \u03c1", title = "Correlation with T / N / stage") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank())

fig5 <- p5a / p5b + plot_layout(heights = c(1.25, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig5_clinical.pdf"), fig5, width = 180, height = 130, units = "mm")
ggsave(file.path(FIG, "Fig5_clinical.png"), fig5, width = 180, height = 130, units = "mm", dpi = 600)

# ---------------- Figure 7: tumor-versus-adjacent discrimination score ----------------
coef <- read.csv(file.path(DATA, "panel_coefficients.csv"))
names(coef) <- c("gene", "coef")
coef <- coef[coef$coef != 0, ]

load_cohort <- function(tag) {
  if (tag == "GSE53625") {
    e <- read.csv(file.path(DATA, "GSE53625_expr.csv"), row.names = 1, check.names = FALSE)
    g <- read.csv(file.path(DATA, "GSE53625_group.csv")); g$sample <- trimws(g$sample)
    e <- e[, g$sample]; y <- as.integer(g$group == "Tumor")
  } else {
    e <- read.csv(file.path(DATA, paste0(tag, "_expr_gene.csv")), row.names = 1, check.names = FALSE)
    g <- read.csv(file.path(DATA, paste0(tag, "_group.csv"))); g$sample <- trimws(g$sample)
    e <- e[, g$sample]; y <- as.integer(g$group == "Tumor")
  }
  e <- as.matrix(e); e <- e[!duplicated(rownames(e)), , drop = FALSE]
  list(expr = e, y = y)
}

# Discovery cohort: nested patient-grouped out-of-fold predictions, NOT resubstitution.
make_roc <- function(y, sc, tag) {
  ord <- order(sc, decreasing = TRUE)
  tp <- cumsum(y[ord] == 1); fp <- cumsum(y[ord] == 0)
  tpr <- c(0, tp / sum(y == 1)); fpr <- c(0, fp / sum(y == 0))
  pos <- sc[y == 1]; neg <- sc[y == 0]
  w <- wilcox.test(pos, neg)$statistic
  auc <- as.numeric(w) / (length(pos) * length(neg))
  list(df = data.frame(fpr = fpr, tpr = tpr,
                       cohort = sprintf("%s (AUC = %.3f)", tag, auc)),
       tab = data.frame(cohort = tag, n_genes = NA_integer_, auc = auc))
}

roc_df <- list(); auc_tab <- list()
oof <- read.csv("D:/escc_sc/data/panel_oof_scores.csv")
r1 <- make_roc(oof$y, oof$oof_score, "GSE53625 nested CV")
roc_df[["GSE53625"]] <- r1$df; auc_tab[["GSE53625"]] <- r1$tab

for (tag in c("GSE20347", "GSE45670")) {
  dd <- load_cohort(tag)
  gs <- intersect(coef$gene, rownames(dd$expr))
  z <- t(scale(t(dd$expr[gs, , drop = FALSE])))
  sc <- as.numeric(crossprod(z, coef$coef[match(gs, coef$gene)]))
  rk <- make_roc(dd$y, sc, tag)
  roc_df[[tag]] <- rk$df; auc_tab[[tag]] <- rk$tab
}
roc <- bind_rows(roc_df); auc_df <- bind_rows(auc_tab)
write.csv(auc_df, file.path(DATA, "panel_auc_R.csv"), row.names = FALSE)

p6a <- ggplot(roc, aes(fpr, tpr, colour = cohort, linetype = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey60") +
  geom_step(linewidth = 0.55) +
  scale_colour_manual(values = c("#b03a2e", "#d98b5f", "#1b6ca8")) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash"), guide = "none") +
  scale_x_continuous(labels = percent, expand = c(0.01, 0)) +
  scale_y_continuous(labels = percent, expand = c(0.01, 0)) +
  labs(x = "1 - specificity", y = "Sensitivity", colour = NULL, title = "Tumor-versus-adjacent discrimination") +
  th + theme(legend.position = c(0.62, 0.18), legend.background = element_blank(),
             legend.key.size = unit(3, "mm"), legend.key.width = unit(8, "mm"),
             legend.text = element_text(size = 6))

c2 <- coef[order(coef$coef), ]
c2$gene <- factor(c2$gene, levels = c2$gene)
p6b <- ggplot(c2, aes(coef, gene, fill = coef > 0)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = c("TRUE" = "#b03a2e", "FALSE" = "#1b6ca8")) +
  labs(x = "Penalized logistic coefficient", y = NULL, title = sprintf("%d-gene descriptive model", nrow(coef))) +
  th + theme(legend.position = "none", axis.text.y = element_text(size = 5.5))

fig6 <- (p6a | p6b) + plot_layout(widths = c(1, 0.85)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig6_panel.pdf"), fig6, width = 175, height = 85, units = "mm")
ggsave(file.path(FIG, "Fig6_panel.png"), fig6, width = 175, height = 85, units = "mm", dpi = 600)
cat("Fig5-6 done\n"); print(auc_df)
