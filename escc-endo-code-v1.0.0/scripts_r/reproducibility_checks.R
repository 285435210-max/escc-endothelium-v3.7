# B4 / M10: clustering-resolution sensitivity, ligand-abundance robustness across
# five subsampling seeds, and a full sessionInfo for the reproducibility record.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix)
})
set.seed(42)
OUT <- "D:/escc_sc/out"
obj <- readRDS(file.path(OUT, "endo_annotated.rds"))

# ---------- 1. clustering resolution 0.4 / 0.6 / 0.8 ----------
base <- as.character(obj$celltype)
ari <- function(x, y) {
  tab <- table(x, y)
  a <- rowSums(tab); b <- colSums(tab); n <- sum(tab)
  c2 <- function(z) sum(z * (z - 1) / 2)
  s1 <- c2(tab); s2 <- c2(a); s3 <- c2(b); s4 <- c2(n)
  (s1 - s2 * s3 / s4) / (0.5 * (s2 + s3) - s2 * s3 / s4)
}
res_rows <- data.frame()
for (r in c(0.4, 0.6, 0.8)) {
  set.seed(42)
  o <- FindClusters(obj, resolution = r, verbose = FALSE)
  lab <- as.character(o$seurat_clusters)
  res_rows <- rbind(res_rows, data.frame(
    resolution = r, n_clusters = length(unique(lab)),
    ARI_vs_9state_annotation = round(ari(lab, base), 4)))
  rm(o); gc(verbose = FALSE)
}
write.csv(res_rows, file.path(OUT, "cluster_resolution_sensitivity.csv"), row.names = FALSE)
print(res_rows)

# ---------- 2. ligand abundance ranking across 5 EC subsampling seeds ----------
obj <- subset(obj, subset = !grepl("^Contam", as.character(celltype)))
genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "ICAM1", "MIF", "APP",
           "CD99", "PECAM1", "CD74", "CXCL9", "CXCL10", "CXCL11",
           "ESM1", "ACKR1", "PLVAP", "SELP", "VCAM1")
genes <- intersect(genes, rownames(obj))
expr <- LayerData(obj, assay = "RNA", layer = "data")
lst <- list()
for (sd in 1:5) {
  set.seed(sd)
  idx <- unlist(lapply(split(colnames(obj), as.character(obj$celltype)),
                       function(x) sample(x, min(length(x), 350))))
  e <- as.matrix(expr[genes, idx, drop = FALSE])
  m <- rowMeans(e, na.rm = TRUE)
  lst[[sd]] <- data.frame(seed = sd, gene = genes, mean_expr = round(m, 3),
                          rank = rank(-m, ties.method = "average"))
}
d <- bind_rows(lst)
write.csv(d, file.path(OUT, "ligand_rank_5seeds.csv"), row.names = FALSE)
sm <- d %>% group_by(gene) %>%
  summarise(mean_rank = mean(rank), min_rank = min(rank), max_rank = max(rank),
            mean_expr = mean(mean_expr), sd_expr = sd(mean_expr), .groups = "drop") %>%
  arrange(mean_rank)
write.csv(sm, file.path(OUT, "ligand_rank_summary.csv"), row.names = FALSE)
print(head(sm, 12))

# ---------- 3. session info ----------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("DONE reproducibility checks\n")
