# Figure 1 & 2: endothelial atlas and compositional remodelling
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})
set.seed(42)
FIG <- "D:/escc_sc/figs"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
OUT <- "D:/escc_sc/out"
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(3, "mm"))

obj <- readRDS(file.path(OUT, "endo_annotated.rds"))
obj <- subset(obj, subset = !grepl("^Contam", as.character(celltype)))
code <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
          "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic")
disp <- c(Arterial = "Arterial", Capillary = "Capillary", Venous = "Venous",
          Venous_inflamed = "Inflamed venule", Tip_angiogenic = "Tip-like",
          Matrix_EC = "Matrix-remodeling", IFNgamma_EC = "Interferon-responsive",
          Proliferating = "Proliferating", Lymphatic = "Lymphatic")
obj$celltype <- factor(obj$celltype, levels = code, labels = unname(disp[code]))
lv <- unname(disp[code])
cols <- c("#1b6ca8","#4a9bc4","#8ec6d9","#c7e2ea","#b03a2e",
          "#d98b5f","#7b4fa0","#4d4d4d","#2e7d32")
names(cols) <- lv

# UMAP labels use short two-line display names to avoid overlap.
short <- c(Arterial = "Arterial", Capillary = "Capillary", Venous = "Venous",
           Venous_inflamed = "Inflamed\nvenule", Tip_angiogenic = "Tip-like",
           Matrix_EC = "Matrix-\nremodeling", IFNgamma_EC = "Interferon-\nresponsive",
           Proliferating = "Proliferating", Lymphatic = "Lymphatic")
back <- setNames(code, unname(disp[code]))
obj$celltype_short <- factor(unname(short[back[as.character(obj$celltype)]]),
                             levels = unname(short[code]))
cols_short <- setNames(unname(cols), unname(short[code]))

# ---- Fig 1a: UMAP ----
p1a <- DimPlot(obj, group.by = "celltype_short", cols = cols_short, pt.size = 0.25, label = TRUE,
               label.size = 2.2, repel = TRUE) +
  labs(title = "ESCC endothelial cells (n = 9,531)", x = "UMAP 1", y = "UMAP 2") +
  th + theme(legend.position = "none")

# ---- Fig 1b: marker dotplot ----
panel <- list(
  Arterial        = c("GJA5","DKK2","HEY2","SEMA3G","FBLN5"),
  Capillary       = c("KIT","APOD","GPIHBP1","RBP7","CA4"),
  Venous          = c("ACKR1","SELP","NR2F2","VCAM1"),
  Venous_inflamed = c("SELE","SELP","VCAM1","ICAM1","SELL","ICAM4","IGFBP5"),
  Tip_angiogenic  = c("ESM1","APLN","NID2","PLVAP"),
  Matrix_EC       = c("IGFBP5","HAPLN1","ADAMTSL2","ITGA8"),
  IFNgamma_EC     = c("CXCL9","CXCL10","CXCL11","RSAD2","IFIT2"),
  Proliferating   = c("MCM10","E2F1","KIF20A","BUB1"),
  Lymphatic       = c("CCL21","PROX1","PDPN","RELN"))
panel <- lapply(panel, function(g) intersect(g, rownames(obj)))
genes <- unique(unlist(panel))
obj$celltype <- factor(obj$celltype, levels = lv)
p1b <- DotPlot(obj, features = genes, group.by = "celltype", dot.scale = 2.6) +
  RotatedAxis() + th +
  scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0) +
  labs(x = NULL, y = NULL, title = "Canonical markers") +
  theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6))

# ---- Fig 1b: QC scatter (all cells, retained vs excluded) ----
suppressPackageStartupMessages(library(Matrix))
C <- readRDS("D:/escc_sc/counts_raw.rds")
mtg <- grepl("^MT-", rownames(C))
nF <- Matrix::colSums(C > 0)
mt <- Matrix::colSums(C[mtg, , drop = FALSE]) / Matrix::colSums(C) * 100
qc <- data.frame(nFeature = nF, mt = mt,
                 status = ifelse(nF >= 200 & nF <= 5000 & mt < 10, "Retained", "Excluded"))
rm(C); gc(verbose = FALSE)
n_ex <- sum(qc$status == "Excluded")
pct_ex <- 100 * n_ex / nrow(qc)

p1c <- ggplot(qc, aes(x = nFeature, y = mt)) +
  geom_point(aes(colour = status), size = 0.22, alpha = 0.55, stroke = 0) +
  geom_vline(xintercept = c(200, 5000), linetype = "dashed", linewidth = 0.25, colour = "grey35") +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.25, colour = "grey35") +
  scale_x_log10(breaks = c(200, 500, 1000, 2000, 5000),
                labels = c("200", "500", "1,000", "2,000", "5,000")) +
  scale_colour_manual(values = c(Retained = "#9db8c9", Excluded = "#b03a2e")) +
  annotate("text", x = 300, y = 19.4, hjust = 0, size = 2.1, colour = "grey20",
           label = sprintf("Before QC  n = %s", format(nrow(qc), big.mark = ","))) +
  annotate("text", x = 300, y = 17.7, hjust = 0, size = 2.1, colour = "#5b7f95",
           label = sprintf("After QC   n = %s", format(sum(qc$status == "Retained"), big.mark = ","))) +
  annotate("text", x = 300, y = 16.0, hjust = 0, size = 2.1, colour = "#b03a2e",
           label = sprintf("Excluded  n = %s (%.1f%%)", format(n_ex, big.mark = ","), pct_ex)) +
  annotate("text", x = 300, y = 14.3, hjust = 0, size = 2.1, colour = "grey20",
           label = "Analyzed*  n = 9,531") +
  labs(x = "Genes per cell", y = "Mitochondrial reads (%)",
       title = "Quality control") +
  th + theme(legend.position = "none")

fig1 <- (p1a | p1c) / p1b + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig1_atlas.pdf"), fig1, width = 180, height = 165, units = "mm")
ggsave(file.path(FIG, "Fig1_atlas.png"), fig1, width = 180, height = 165, units = "mm", dpi = 600)

# ---- Fig 2: composition ----
prop <- read.csv("D:/escc_sc/data/celltype_proportion_closed.csv")
prop <- prop %>% filter(!grepl("^Contam", celltype)) %>%
  mutate(celltype = factor(celltype, levels = code, labels = lv),
         patient = sub("T$", "", sample))
# 未检出的状态记为 0，构成完整样本×状态网格
all_s <- unique(prop[, c("sample", "type")])
grid <- expand.grid(sample = all_s$sample, celltype = lv, stringsAsFactors = FALSE)
grid$type <- all_s$type[match(grid$sample, all_s$sample)]
prop <- grid %>%
  left_join(prop[, c("sample", "type", "celltype", "n", "frac")],
            by = c("sample", "type", "celltype")) %>%
  mutate(frac = ifelse(is.na(frac), 0, frac), n = ifelse(is.na(n), 0, n),
         patient = sub("T$", "", sample),
         celltype = factor(as.character(celltype), levels = lv))

# 按静脉型占比降序排列，使梯度可见
ord_s <- prop %>% filter(type == "Tumor", celltype == "Venous") %>%
  arrange(desc(frac)) %>% pull(sample)
p2a <- prop %>% filter(type == "Tumor") %>%
  mutate(sample = factor(sample, levels = ord_s)) %>%
  ggplot(aes(x = sample, y = frac, fill = celltype)) +
  geom_col(width = 1) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(x = "Tumor samples (n = 60)", y = "Fraction of endothelial cells", fill = NULL,
       title = "Per-patient endothelial composition") +
  th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
             legend.position = "right", legend.text = element_text(size = 6))

# Four tumors are paired with the four normal libraries; the primary
# independent comparison uses the remaining 56 tumors (Table S12).
paired_tumors <- c("P126T", "P127T", "P128T", "P130T")
key <- lv
p2b <- prop %>% filter(celltype %in% key,
                        !(type == "Tumor" & sample %in% paired_tumors)) %>%
  ggplot(aes(x = celltype, y = frac, fill = type)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.62,
               position = position_dodge(0.75)) +
  geom_point(aes(colour = type), size = 0.55, alpha = 0.55, stroke = 0,
             position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.75)) +
  scale_colour_manual(values = c(Normal = "#2f6d8f", Tumor = "#7e2418"), guide = "none") +
  scale_fill_manual(values = c(Normal = "#4a9bc4", Tumor = "#b03a2e"),
                     labels = c(Normal = "Adjacent non-malignant", Tumor = "Tumor")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Fraction", fill = NULL,
       title = "Unmatched tumors vs adjacent non-malignant mucosa") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# median fraction comparison (tumour vs normal)
med <- prop %>% filter(celltype %in% key,
                       !(type == "Tumor" & sample %in% paired_tumors)) %>%
  group_by(celltype, type) %>% summarise(m = median(frac), .groups = "drop")
p2c <- ggplot(med, aes(x = celltype, y = m, fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.62, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Normal = "#4a9bc4", Tumor = "#b03a2e"),
                     labels = c(Normal = "Adjacent non-malignant", Tumor = "Tumor")) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Median fraction", fill = NULL, title = "Median composition shift") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

fig2 <- (p2a / (p2b | p2c)) + plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig2_composition.pdf"), fig2, width = 200, height = 160, units = "mm")
ggsave(file.path(FIG, "Fig2_composition.png"), fig2, width = 200, height = 160, units = "mm", dpi = 600)

cat("Fig1-2 done\n")
