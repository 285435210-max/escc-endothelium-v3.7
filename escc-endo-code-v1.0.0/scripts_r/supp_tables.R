suppressPackageStartupMessages({library(Seurat); library(dplyr); library(Matrix)})
OUT <- "D:/escc_sc/supp"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS("D:/escc_sc/out/endo_annotated.rds")
obj <- subset(obj, subset = !grepl("^Contam", as.character(celltype)))
obj$celltype <- factor(as.character(obj$celltype))
Idents(obj) <- "celltype"
cat("cells:", ncol(obj), "\n")

MKF <- file.path(OUT, "TableS2_state_markers.csv")
if (file.exists(MKF)) { cat("markers cached\n") } else {
mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
mk <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 40) %>% arrange(cluster, desc(avg_log2FC))
write.csv(mk, file.path(OUT, "TableS2_state_markers.csv"), row.names = FALSE)
cat("markers:", nrow(mk), "\n")
}

# ---- QC per library ----
raw <- readRDS("D:/escc_sc/counts_raw.rds")
hdr <- readLines(gzfile("D:/龙虾文档/2026-11-02-09-16-41/ESCC_endo/data_manifest/GSE160269_UMI_matrix_Endothelial.txt.gz"), n = 1)
bc <- strsplit(hdr, " ")[[1]]; bc <- bc[bc != ""]
colnames(raw) <- bc
mtg <- grepl("^MT-", rownames(raw))
nf <- Matrix::colSums(raw > 0)
mt <- Matrix::colSums(raw[mtg, , drop = FALSE]) / Matrix::colSums(raw) * 100
smp <- sub("-[^-]+$", "", colnames(raw))
pre <- data.frame(sample = smp, nFeature = nf, mt = mt) %>%
  mutate(pass = nFeature >= 200 & nFeature <= 5000 & mt < 10) %>%
  group_by(sample) %>%
  summarise(n_cells_raw = n(), n_cells_qc_pass = sum(pass),
            median_genes_raw = median(nFeature),
            median_mt_raw = round(median(mt), 2), .groups = "drop")
post <- data.frame(sample = sub("-[^-]+$", "", colnames(obj)),
                   nFeature = obj$nFeature_RNA, mt = obj$percent.mt) %>%
  group_by(sample) %>%
  summarise(n_cells_qc = n(), median_genes_qc = median(nFeature),
            median_mt_qc = round(median(mt), 2), .groups = "drop")
qc <- full_join(pre, post, by = "sample") %>%
  mutate(n_cells_qc = ifelse(is.na(n_cells_qc), 0, n_cells_qc),
         pct_retained = round(100 * n_cells_qc / n_cells_raw, 1)) %>%
  arrange(sample)
write.csv(qc, file.path(OUT, "TableS1_QC_by_library.csv"), row.names = FALSE)
cat("libraries:", nrow(qc), " total raw:", sum(qc$n_cells_raw), " qc:", sum(qc$n_cells_qc), "\n")
