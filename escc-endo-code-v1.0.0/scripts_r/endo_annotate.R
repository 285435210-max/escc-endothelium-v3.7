# Annotate endothelial clusters and quantify populations
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "D:/escc_sc/out"
obj <- readRDS(file.path(OUT, "endo_clustered.rds"))

map <- c(
  "0"  = "Venous",
  "1"  = "Tip_angiogenic",
  "2"  = "Matrix_EC",
  "3"  = "Venous_inflamed",
  "4"  = "Venous",
  "5"  = "Arterial",
  "6"  = "Venous",
  "7"  = "Capillary",
  "8"  = "Contam_fibroblast",
  "9"  = "Arterial",
  "10" = "IFNgamma_EC",
  "11" = "Proliferating",
  "12" = "Proliferating",
  "13" = "Capillary",
  "14" = "Capillary",
  "15" = "Contam_pericyte",
  "16" = "Contam_immune",
  "17" = "Lymphatic"
)
ct <- factor(unname(map[as.character(obj$seurat_clusters)]),
             levels = c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
                        "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic",
                        "Contam_pericyte","Contam_fibroblast","Contam_immune"))
names(ct) <- colnames(obj)
obj$celltype <- ct
cat("celltype counts:\n"); print(table(obj$celltype))

png(file.path(OUT, "UMAP_celltype.png"), width = 2000, height = 1500, res = 300)
print(DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.3) +
        ggtitle("Endothelial populations in ESCC (GSE160269)"))
dev.off()

png(file.path(OUT, "UMAP_celltype_split.png"), width = 2400, height = 1000, res = 300)
print(DimPlot(obj, group.by = "celltype", split.by = "type", label = TRUE, repel = TRUE,
              pt.size = 0.3, ncol = 2))
dev.off()

# ---------- proportions per sample ----------
tab <- obj@meta.data %>% count(sample, type, celltype)
tot <- tab %>% group_by(sample) %>% summarise(n_total = sum(n), .groups = "drop")
prop <- tab %>% left_join(tot, by = "sample") %>% mutate(frac = n / n_total)
write.csv(prop, file.path(OUT, "celltype_proportion_by_sample.csv"), row.names = FALSE)

summ <- prop %>% group_by(celltype, type) %>%
  summarise(median_frac = median(frac), mean_frac = mean(frac), n_samples = n(), .groups = "drop")
write.csv(summ, file.path(OUT, "celltype_proportion_summary.csv"), row.names = FALSE)
cat("\nmedian fraction by type:\n")
print(as.data.frame(summ %>% select(celltype, type, median_frac) %>%
                      pivot_wider(names_from = type, values_from = median_frac)))

png(file.path(OUT, "Proportion_by_type.png"), width = 2200, height = 1100, res = 300)
print(ggplot(prop, aes(x = celltype, y = frac, fill = type)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
        scale_fill_manual(values = c(Normal = "#2e6da4", Tumor = "#b03a2e")) +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = NULL, y = "Fraction of endothelial cells per sample", fill = NULL))
dev.off()

# ---------- tumor vs normal per patient (tumor-dominant; normal n=4) ----------
cells_keep <- c("Tip_angiogenic", "IFNgamma_EC", "Proliferating", "Lymphatic", "Arterial", "Capillary", "Venous")
p <- prop %>% filter(celltype %in% cells_keep) %>%
  group_by(celltype) %>%
  summarise(
    n_tumor = sum(type == "Tumor"), n_normal = sum(type == "Normal"),
    med_tumor = median(frac[type == "Tumor"]),
    med_normal = median(frac[type == "Normal"]),
    p_wilcox = tryCatch(wilcox.test(frac[type == "Tumor"], frac[type == "Normal"])$p.value, error = function(e) NA),
    .groups = "drop")
write.csv(p, file.path(OUT, "celltype_tumor_vs_normal.csv"), row.names = FALSE)
cat("\ntumor vs normal (patient-level fractions):\n"); print(as.data.frame(p))

saveRDS(obj, file.path(OUT, "endo_annotated.rds"))
cat("\nDONE annotation\n")
