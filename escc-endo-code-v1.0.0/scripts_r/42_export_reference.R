# Export GSE160269 endothelial reference (9 states) for external label transfer
# NOTE: all output paths ASCII-only (Windows GBK locale cannot write CJK paths from R)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
})
REF <- "D:/escc_sc/out/endo_annotated.rds"
OUT <- "D:/escc_sc/external"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

STATES <- c("Arterial", "Capillary", "Venous", "Venous_inflamed",
            "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC",
            "Proliferating", "Lymphatic")

obj <- readRDS(REF)
cat("ref dim:", dim(obj), "\n")
lab <- as.character(obj@meta.data$celltype)
names(lab) <- colnames(obj)
print(table(lab))

keep <- lab %in% STATES
cat("cells kept:", sum(keep), "/", length(lab), "\n")

ex <- GetAssayData(obj, layer = "data")[, keep, drop = FALSE]
ex <- as(ex, "dgCMatrix")
lab <- lab[colnames(ex)]

# 1. labels
write.csv(data.frame(cell = colnames(ex), label = lab,
                     sample = as.character(obj@meta.data[colnames(ex), "sample"]),
                     type = as.character(obj@meta.data[colnames(ex), "type"])),
          file.path(OUT, "ref_labels.csv"), row.names = FALSE)

# 2. sparse matrix (genes x cells) in MatrixMarket format
Matrix::writeMM(ex, file.path(OUT, "ref_matrix.mtx"))
writeLines(rownames(ex), file.path(OUT, "ref_genes.txt"))
writeLines(colnames(ex), file.path(OUT, "ref_barcodes.txt"))
cat("reference expr dim:", dim(ex), " nonzero:", length(ex@x), "\n")

# 3. state centroid (mean log-normalised expression per state)
cent <- do.call(cbind, lapply(STATES, function(s) {
  Matrix::rowMeans(ex[, lab == s, drop = FALSE])
}))
colnames(cent) <- STATES
write.csv(data.frame(gene = rownames(ex), cent, check.names = FALSE),
          file.path(OUT, "ref_state_centroids.csv"), row.names = FALSE)

# 4. top-20 positive markers per state (module-score fallback)
m <- fread("D:/escc_sc/supp/TableS2_state_markers.csv")
m <- m[cluster %in% STATES & avg_log2FC > 0]
m <- m[order(cluster, -avg_log2FC)]
top <- m[, .SD[1:min(20, .N)], by = cluster]
write.csv(top[, .(cluster, gene, avg_log2FC, p_val_adj)],
          file.path(OUT, "ref_state_markers_top20.csv"), row.names = FALSE)
cat("marker rows:", nrow(top), " states:", length(unique(top$cluster)), "\n")
cat("\nDONE\n")
