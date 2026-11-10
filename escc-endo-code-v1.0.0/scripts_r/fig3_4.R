# Figure 3 & 4: tip-cell program and interferon-responsive endothelium
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggrepel)
})
FIG <- "D:/escc_sc/figs"; DATA <- "D:/escc_sc/data"; OUT <- "D:/escc_sc/out"
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 9, face = "bold"))

obj <- readRDS(file.path(OUT, "endo_annotated.rds"))
obj <- subset(obj, subset = !grepl("^Contam", as.character(celltype)))
lv <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
        "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic")
disp <- c(Arterial = "Arterial", Capillary = "Capillary", Venous = "Venous",
          Venous_inflamed = "Inflamed venule", Tip_angiogenic = "Tip-like",
          Matrix_EC = "Matrix-remodeling", IFNgamma_EC = "Interferon-responsive",
          Proliferating = "Proliferating", Lymphatic = "Lymphatic")
obj$celltype <- factor(obj$celltype, levels = lv, labels = unname(disp[lv]))

# ---------------- volcano: Tip vs Venous ----------------
de <- read.csv(file.path(DATA, "DE_Tip_vs_Venous.csv"))
de$sig <- ifelse(de$adj.P.Val < 0.05 & de$logFC > 1, "Up in tip",
          ifelse(de$adj.P.Val < 0.05 & de$logFC < -1, "Down in tip", "NS"))
lab_up <- de %>% filter(sig == "Up in tip") %>% arrange(desc(logFC)) %>% head(9)
lab_dn <- de %>% filter(sig == "Down in tip") %>% arrange(logFC) %>% head(5)
p3a <- ggplot(de, aes(logFC, -log10(P.Value), colour = sig)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_colour_manual(values = c("Up in tip" = "#b03a2e", "Down in tip" = "#1b6ca8", "NS" = "grey85")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) +
  geom_text_repel(data = lab_up, aes(label = gene), size = 1.9, max.overlaps = 30,
                  segment.size = 0.12, min.segment.length = 0, nudge_x = 1.2,
                  direction = "y", hjust = 0) +
  geom_text_repel(data = lab_dn, aes(label = gene), size = 1.9, max.overlaps = 30,
                  segment.size = 0.12, min.segment.length = 0, nudge_x = -1.2,
                  direction = "y", hjust = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = expression(log[2]~"fold change (tip vs venous)"), y = expression(-log[10]~italic(P)),
       colour = NULL, title = "Tip cells vs venous endothelium") +
  th + theme(legend.position = "top", legend.direction = "horizontal",
             legend.key.size = unit(2.5, "mm"))

# ---------------- enrichGO: tip signature ----------------
tip <- read.csv(file.path(DATA, "Tip_signature_shared.csv"))$gene
eg <- bitr(tip, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_tip <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_tip_df <- as.data.frame(go_tip)
write.csv(go_tip_df, file.path(DATA, "enrichGO_tip.csv"), row.names = FALSE)
cat("tip GO terms:", nrow(go_tip_df), "\n")
gd <- as.data.frame(go_tip) %>% slice_head(n = 12) %>%
  mutate(gr = vapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), 0),
         nlp = -log10(p.adjust),
         Description = factor(Description, levels = rev(Description)))
p3b <- ggplot(gd, aes(gr, Description)) +
  geom_point(aes(size = Count, colour = nlp)) +
  scale_colour_gradient(low = "#9ec4d8", high = "#b03a2e",
                        name = expression(-log[10]~italic(P)[adj])) +
  scale_size_continuous(name = "Count", range = c(0.8, 4.2)) +
  labs(x = "Gene ratio", y = NULL, title = "GO biological processes (tip signature)") +
  th + theme(axis.text.y = element_text(size = 6),
             legend.key.size = unit(2.6, "mm"),
             legend.text = element_text(size = 6),
             legend.title = element_text(size = 6.5))

# ---------------- Fig 4a: interferon-responsive EC markers ----------------
mk <- c("CXCL9","CXCL10","CXCL11","RSAD2","IFIT2","HLA-DRA","ICAM1","SELE","VWF","ACKR1")
mk <- intersect(mk, rownames(obj))
p4a <- DotPlot(obj, features = mk, group.by = "celltype", dot.scale = 2.8) + RotatedAxis() +
  scale_colour_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0) +
  labs(x = NULL, y = NULL, title = "Inflammatory / antigen-presenting program") + th +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right", legend.box = "vertical",
        legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm"),
        legend.text = element_text(size = 5.5), legend.title = element_text(size = 6))

# ---------------- Fig 4b: enrichGO of interferon-responsive up genes ----------------
ifn <- read.csv(file.path(DATA, "DE_IFNg_vs_Venous.csv"))
ifn_up <- ifn$gene[ifn$adj.P.Val < 0.05 & ifn$logFC > 1]
eg2 <- bitr(ifn_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_ifn <- enrichGO(eg2$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_ifn_df <- as.data.frame(go_ifn)
write.csv(go_ifn_df, file.path(DATA, "enrichGO_IFNg.csv"), row.names = FALSE)
cat("IFNg GO terms:", nrow(go_ifn_df), "\n")
gi <- as.data.frame(go_ifn) %>% slice_head(n = 12) %>%
  mutate(gr = vapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), 0),
         nlp = -log10(p.adjust),
         Description = factor(Description, levels = rev(Description)))
p4b <- ggplot(gi, aes(gr, Description)) +
  geom_point(aes(size = Count, colour = nlp)) +
  scale_colour_gradient(low = "#9ec4d8", high = "#b03a2e",
                        name = expression(-log[10]~italic(P)[adj])) +
  scale_size_continuous(name = "Count", range = c(0.8, 4.2)) +
  labs(x = "Gene ratio", y = NULL, title = "GO biological processes") +
  th + theme(axis.text.y = element_text(size = 6),
             legend.key.size = unit(2.6, "mm"),
             legend.text = element_text(size = 6),
             legend.title = element_text(size = 6.5))

# ---------------- Fig 4c: interferon-responsive fraction per tumour ----------------
# NB: celltype_proportion_closed.csv stores machine-code state names.
prop <- read.csv(file.path(DATA, "celltype_proportion_closed.csv")) %>%
  filter(type == "Tumor") %>% mutate(patient = sub("T$", "", sample))
d <- data.frame(patient = sort(unique(prop$patient))) %>%
  left_join(prop %>% filter(celltype == "IFNgamma_EC") %>% transmute(patient, frac),
            by = "patient") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac))
p4c <- ggplot(d, aes(x = reorder(patient, -frac), y = frac)) +
  geom_col(fill = "#7b4fa0", width = 0.85) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Tumor samples", y = "Interferon-responsive fraction",
       title = "Prevalence across tumors") + th +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig3 <- p3a / p3b + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig3_tip.pdf"), fig3, width = 170, height = 175, units = "mm")
ggsave(file.path(FIG, "Fig3_tip.png"), fig3, width = 170, height = 175, units = "mm", dpi = 600)

fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1, 0.62)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig4_ifng.pdf"), fig4, width = 180, height = 175, units = "mm")
ggsave(file.path(FIG, "Fig4_ifng.png"), fig4, width = 180, height = 175, units = "mm", dpi = 600)
cat("Fig3-4 done\n")
