# B4 robustness, v3.5: nboot=200, two EC subsampling seeds.
# T-cell labels now include an explicit "Other/unclassified T" group.
# Outputs go to D:/escc_sc/cellchat_out/robust_v35 so v3.4 results are preserved.
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(Matrix); library(dplyr)
})
options(future.globals.maxSize = 4 * 1024^3)
OUT <- "D:/escc_sc/cellchat_out/robust_v35"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(OUT, "run.log")
cat("START", format(Sys.time()), "\n", file = logf, append = FALSE)

ec_all <- readRDS("D:/escc_sc/out/endo_annotated.rds")
ec_all <- subset(ec_all, subset = !grepl("^Contam", as.character(celltype)))
ec_all$grp <- as.character(ec_all$celltype)

tc <- readMM("D:/escc_sc/cellchat/tcell_subset.mtx")
tg <- readLines("D:/escc_sc/cellchat/tcell_genes.txt")
tb <- readLines("D:/escc_sc/cellchat/tcell_barcodes.txt")
rownames(tc) <- tg; colnames(tc) <- tb
tc <- as(tc, "CsparseMatrix")
det <- function(g) if (g %in% rownames(tc)) tc[g, ] > 0 else rep(FALSE, ncol(tc))
cd8 <- det("CD8A") | det("CD8B")
cd4 <- det("CD4") & !cd8
treg <- cd4 & det("FOXP3") & det("IL2RA")
exf <- det("PDCD1") + det("HAVCR2") + det("LAG3") + det("TIGIT") + det("TOX")
tex <- exf >= 2
nk <- (!det("CD3D") & !det("CD3E")) & (det("GNLY") | det("NKG7") | det("KLRD1"))
lab <- rep("Other/unclassified T", ncol(tc))
lab[cd4 & !treg] <- "CD4 T"
lab[cd8] <- "CD8 T"
lab[cd8 & tex] <- "CD8 exhausted"
lab[treg] <- "Treg"
lab[nk] <- "NK-like"
names(lab) <- colnames(tc)
write.csv(as.data.frame(table(lab)), file.path(OUT, "tcell_group_counts.csv"),
          row.names = FALSE)

ec9 <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
         "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic")
tt <- c("CD8 T","CD8 exhausted","CD4 T","Treg","NK-like")

for (sd in c(1, 2)) {
  cat("seed", sd, "start", format(Sys.time()), "\n", file = logf, append = TRUE)
  set.seed(sd)
  idx <- unlist(lapply(split(colnames(ec_all), ec_all$grp),
                       function(x) sample(x, min(length(x), 350))))
  ec <- subset(ec_all, cells = idx)
  ec_counts <- LayerData(ec, assay = "RNA", layer = "counts")
  common <- intersect(rownames(ec_counts), rownames(tc))
  mat <- cbind(ec_counts[common, ], tc[common, ])
  grp <- c(ec$grp, lab[colnames(tc)])
  meta <- data.frame(labels = grp,
                     samples = c(as.character(ec$sample), sub("-.*$", "", colnames(tc))),
                     row.names = colnames(mat))
  obj <- CreateSeuratObject(counts = mat, meta.data = meta)
  obj <- NormalizeData(obj, verbose = FALSE)
  data.input <- LayerData(obj, assay = "RNA", layer = "data")
  rm(mat, obj); gc(verbose = FALSE)

  cc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", nboot = 200, seed.use = sd)
  cc <- filterCommunication(cc, min.cells = 10)
  df <- subsetCommunication(cc)
  df$fdr <- p.adjust(df$pval, method = "BH")
  write.csv(df, file.path(OUT, sprintf("net_seed%d_nboot200.csv", sd)), row.names = FALSE)

  s <- df[df$source %in% ec9 & df$target %in% tt, ]
  pw <- aggregate(prob ~ pathway_name, s, sum)
  pw <- pw[order(-pw$prob), ]
  write.csv(pw, file.path(OUT, sprintf("pathway_ECtoT_seed%d.csv", sd)), row.names = FALSE)

  smry <- data.frame(
    seed = sd, nboot = 200, n_interactions = nrow(df),
    pval_zero = sum(df$pval == 0),
    min_nonzero_pval = if (any(df$pval > 0)) min(df$pval[df$pval > 0]) else NA,
    nominal_sig = sum(df$pval < 0.05), fdr_sig = sum(df$fdr < 0.05),
    n_ECtoT = nrow(s),
    top1 = pw$pathway_name[1], top2 = pw$pathway_name[2], top3 = pw$pathway_name[3])
  write.csv(smry, file.path(OUT, sprintf("summary_seed%d.csv", sd)), row.names = FALSE)
  cat("DONE seed", sd, format(Sys.time()), "\n", file = logf, append = TRUE)
  rm(cc, df, data.input, ec, ec_counts, pw, s); gc(verbose = FALSE)
}
cat("ALL DONE", format(Sys.time()), "\n", file = logf, append = TRUE)
