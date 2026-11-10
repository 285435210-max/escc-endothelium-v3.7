# v3.7 reviewer question 1: is cluster 3 (Venous_inflamed) genuinely endothelial?
# Quantifies pan-endothelial, activation, leukocyte and erythroid markers per state
# and per cell within the Venous_inflamed population.
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})
OUT <- "D:/escc_sc/out_v37"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS("D:/escc_sc/out/endo_annotated.rds")
mat <- GetAssayData(obj, assay = "RNA", layer = "data")
ct <- as.character(obj$celltype)

PANEL <- c("PECAM1", "CDH5", "VWF", "KDR", "EMCN", "CLDN5", "EGFL7",
           "SELL", "SELP", "SELE", "ICAM4", "ICAM1", "IGFBP5",
           "ACKR1", "VCAM1", "NR2F2", "APLNR", "ESM1", "CXCL9",
           "PTPRC", "LST1", "TYROBP", "CD3D", "CD8A", "LYZ",
           "HBA1", "HBA2", "HBB", "GYPA", "ALAS2", "EPB42")
g <- PANEL[PANEL %in% rownames(mat)]
groups <- c("Arterial", "Capillary", "Venous", "Venous_inflamed",
            "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC", "Proliferating",
            "Lymphatic", "Contam_pericyte", "Contam_fibroblast", "Contam_immune")

res <- do.call(rbind, lapply(groups, function(s) {
  idx <- which(ct == s)
  if (!length(idx)) return(NULL)
  sub <- as.matrix(mat[g, idx, drop = FALSE])
  data.frame(state = s, n_cells = length(idx), gene = g,
             pct_detected = round(100 * rowSums(sub > 0) / ncol(sub), 1),
             mean_lognorm = round(rowMeans(sub), 3),
             stringsAsFactors = FALSE)
}))
write.csv(res, file.path(OUT, "inflamed_marker_by_state.csv"), row.names = FALSE)

PAN <- intersect(c("PECAM1", "CDH5", "VWF", "KDR", "EMCN", "CLDN5", "EGFL7"), rownames(mat))
ACT <- intersect(c("SELE", "SELP", "VCAM1", "NR2F2", "APLNR", "ICAM1"), rownames(mat))
LEU <- intersect(c("PTPRC", "LST1", "TYROBP", "CD3D", "CD8A", "LYZ"), rownames(mat))
ERY <- intersect(c("HBA1", "HBA2", "HBB", "GYPA", "ALAS2", "EPB42"), rownames(mat))

per_state <- do.call(rbind, lapply(groups, function(s) {
  idx <- which(ct == s)
  if (!length(idx)) return(NULL)
  n_pan <- colSums(mat[PAN, idx, drop = FALSE] > 0)
  n_act <- colSums(mat[ACT, idx, drop = FALSE] > 0)
  n_leu <- colSums(mat[LEU, idx, drop = FALSE] > 0)
  n_ery <- colSums(mat[ERY, idx, drop = FALSE] > 0)
  data.frame(state = s, n_cells = length(idx),
             pct_with_ge1_panEC = round(100 * mean(n_pan >= 1), 1),
             pct_with_ge2_panEC = round(100 * mean(n_pan >= 2), 1),
             median_n_panEC = median(n_pan),
             pct_with_ge2_activation = round(100 * mean(n_act >= 2), 1),
             pct_with_ge1_leukocyte = round(100 * mean(n_leu >= 1), 1),
             pct_with_ge1_erythroid = round(100 * mean(n_ery >= 1), 1),
             stringsAsFactors = FALSE)
}))
write.csv(per_state, file.path(OUT, "inflamed_per_cell_check.csv"), row.names = FALSE)

idx <- which(ct == "Venous_inflamed")
n_pan <- colSums(mat[PAN, idx, drop = FALSE] > 0)
n_act <- colSums(mat[ACT, idx, drop = FALSE] > 0)
n_leu <- colSums(mat[LEU, idx, drop = FALSE] > 0)
n_ery <- colSums(mat[ERY, idx, drop = FALSE] > 0)
sell <- if ("SELL" %in% rownames(mat)) as.numeric(mat["SELL", idx] > 0) else rep(NA_real_, length(idx))
icam4 <- if ("ICAM4" %in% rownames(mat)) as.numeric(mat["ICAM4", idx] > 0) else rep(NA_real_, length(idx))
co <- data.frame(
  subset = c("all Venous_inflamed", "SELL+", "ICAM4+", "SELL+ or ICAM4+",
             ">=2 pan-endothelial", ">=3 pan-endothelial"),
  n_cells = c(length(idx), sum(sell == 1, na.rm = TRUE), sum(icam4 == 1, na.rm = TRUE),
              sum((sell == 1 | icam4 == 1), na.rm = TRUE), sum(n_pan >= 2), sum(n_pan >= 3)),
  pct_with_ge2_panEC = round(c(100 * mean(n_pan >= 2),
    100 * mean(n_pan[sell == 1] >= 2, na.rm = TRUE),
    100 * mean(n_pan[icam4 == 1] >= 2, na.rm = TRUE),
    100 * mean(n_pan[(sell == 1 | icam4 == 1)] >= 2, na.rm = TRUE),
    100, 100), 1),
  pct_with_ge2_activation = round(c(100 * mean(n_act >= 2),
    100 * mean(n_act[sell == 1] >= 2, na.rm = TRUE),
    100 * mean(n_act[icam4 == 1] >= 2, na.rm = TRUE),
    100 * mean(n_act[(sell == 1 | icam4 == 1)] >= 2, na.rm = TRUE),
    100 * mean(n_act[n_pan >= 2] >= 2),
    100 * mean(n_act[n_pan >= 3] >= 2)), 1),
  pct_with_ge1_leukocyte = round(c(100 * mean(n_leu >= 1),
    100 * mean(n_leu[sell == 1] >= 1, na.rm = TRUE),
    100 * mean(n_leu[icam4 == 1] >= 1, na.rm = TRUE),
    100 * mean(n_leu[(sell == 1 | icam4 == 1)] >= 1, na.rm = TRUE),
    100 * mean(n_leu[n_pan >= 2] >= 1),
    100 * mean(n_leu[n_pan >= 3] >= 1)), 1),
  pct_with_ge1_erythroid = round(c(100 * mean(n_ery >= 1),
    100 * mean(n_ery[sell == 1] >= 1, na.rm = TRUE),
    100 * mean(n_ery[icam4 == 1] >= 1, na.rm = TRUE),
    100 * mean(n_ery[(sell == 1 | icam4 == 1)] >= 1, na.rm = TRUE),
    100 * mean(n_ery[n_pan >= 2] >= 1),
    100 * mean(n_ery[n_pan >= 3] >= 1)), 1))
write.csv(co, file.path(OUT, "inflamed_coexpression.csv"), row.names = FALSE)

cat("=== per-state pan-endothelial / lineage check ===\n")
print(per_state, row.names = FALSE)
cat("\n=== Venous_inflamed co-expression ===\n")
print(co, row.names = FALSE)
cat("\n=== Venous_inflamed marker detection (%) ===\n")
print(res[res$state == "Venous_inflamed", c("gene", "pct_detected", "mean_lognorm")], row.names = FALSE)
cat("\nDONE identity check\n")
