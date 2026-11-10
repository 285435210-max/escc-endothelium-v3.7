# GSE160269 endothelial cells: QC, harmony, clustering, annotation
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat)
  library(ggplot2); library(dplyr); library(harmony)
})

TSV <- "D:/escc_sc/endo.tsv"
RDS <- "D:/escc_sc/counts_raw.rds"
OUT <- "D:/escc_sc/out"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

con <- file(TSV, "rt"); hdr <- readLines(con, n = 1); close(con)
cells <- strsplit(hdr, " ")[[1]]; cells <- cells[cells != ""]
cat("cells:", length(cells), "\n")

if (file.exists(RDS)) {
  counts <- readRDS(RDS)
  cat("loaded cached counts\n")
} else {
  chunk <- 1500; start <- 2L; blocks <- list()
  repeat {
    dt <- tryCatch(fread(TSV, sep = " ", header = FALSE, skip = start - 1L,
                         nrows = chunk, showProgress = FALSE),
                   error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) break
    g <- dt[[1]]; m <- as.matrix(dt[, -1, with = FALSE]); storage.mode(m) <- "double"
    keep <- Matrix::rowSums(m > 0) >= 3
    if (any(keep)) {
      sm <- Matrix(m[keep, , drop = FALSE], sparse = TRUE); rownames(sm) <- g[keep]
      blocks[[length(blocks) + 1L]] <- sm
    }
    start <- start + nrow(dt); rm(dt, m); gc(verbose = FALSE)
  }
  counts <- do.call(rbind, blocks); rm(blocks); gc(verbose = FALSE)
  saveRDS(counts, RDS)
}
colnames(counts) <- cells
cat("counts dim:", dim(counts), " nonzeros:", length(counts@x), "\n")

meta <- data.frame(cell = cells, stringsAsFactors = FALSE)
meta$sample <- sub("^(P[0-9]+[TN]).*$", "\\1", meta$cell)
meta$type   <- ifelse(grepl("T$", meta$sample), "Tumor", "Normal")
rownames(meta) <- meta$cell
cat("samples:", length(unique(meta$sample)), " tumor:", sum(meta$type == "Tumor"),
    " normal:", sum(meta$type == "Normal"), "\n")

QC_RDS <- file.path(OUT, "endo_qc.rds")
if (file.exists(QC_RDS)) {
  obj <- readRDS(QC_RDS); cat("loaded cached QC object:", ncol(obj), "cells\n")
} else {
obj <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 0, min.features = 200)
rm(counts); gc(verbose = FALSE)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
cat("after min.features:", ncol(obj), "cells\n")

obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt < 10)
cat("after QC:", ncol(obj), "cells\n")
print(summary(obj$nFeature_RNA))
saveRDS(obj, file.path(OUT, "endo_qc.rds"))
}

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
saveRDS(obj, file.path(OUT, "endo_clustered.rds"))
cat("clusters:\n"); print(table(obj$seurat_clusters))

mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(mk, file.path(OUT, "cluster_markers.csv"), row.names = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% ungroup()
write.csv(top, file.path(OUT, "cluster_top5_markers.csv"), row.names = FALSE)

panel <- list(
  Pan_endothelial = c("PECAM1","CDH5","VWF","CLDN5","EGFL7"),
  Arterial        = c("GJA5","SEMA3G","HEY1","DKK2","LTBP4"),
  Venous          = c("ACKR1","SELP","VCAM1","NR2F2"),
  Capillary       = c("CA4","RGCC","CD36","MMRN2","BTNL9"),
  Tip_angiogenic  = c("ESM1","PLVAP","APLN","COL4A1","INSR"),
  Lymphatic       = c("LYVE1","CCL21","PROX1","PDPN","MMP1"),
  Pericyte_contam = c("RGS5","PDGFRB","ACTA2","MCAM","NOTCH3")
)
panel <- lapply(panel, function(g) intersect(g, rownames(obj)))

png(file.path(OUT, "UMAP_clusters.png"), width = 1800, height = 1400, res = 300)
print(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
        ggtitle("GSE160269 endothelial cells"))
dev.off()

png(file.path(OUT, "UMAP_type.png"), width = 2200, height = 900, res = 300)
print(DimPlot(obj, reduction = "umap", group.by = "type", split.by = "type",
              cols = c("#b03a2e", "#2e6da4")))
dev.off()

png(file.path(OUT, "DotPlot_panel.png"), width = 2200, height = 1200, res = 300)
print(DotPlot(obj, features = panel, group.by = "seurat_clusters") + RotatedAxis())
dev.off()

comp <- as.data.frame(table(obj$seurat_clusters, obj$type))
colnames(comp) <- c("cluster", "type", "n")
write.csv(comp, file.path(OUT, "cluster_composition.csv"), row.names = FALSE)

allg <- unique(unlist(panel))
av <- as.matrix(AverageExpression(obj, features = allg, group.by = "seurat_clusters", slot = "data")$RNA)
write.csv(round(av, 3), file.path(OUT, "panel_mean_by_cluster.csv"))
cat("\n=== panel mean by cluster ===\n")
print(round(av, 2))

cat("\nDONE\n")
