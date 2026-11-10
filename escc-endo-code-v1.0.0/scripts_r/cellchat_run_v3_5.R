# CellChat v2: endothelial - T cell communication (GSE160269)
# v3.5 change: marker-negative T cells are NO LONGER assigned to "CD4 T" by default.
# They receive an explicit "Other/unclassified T" label. The five interpretable
# T-cell groups are unchanged in definition; CD8 exhausted requires >=2 markers.
# All downstream EC->T summaries exclude the unclassified group, and its size is
# reported in tcell_group_counts_v3_5.csv.
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(Matrix); library(ggplot2); library(dplyr)
})
set.seed(42)
OUT <- "D:/escc_sc/cellchat_out"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 4 * 1024^3)

# ---------- endothelial (subsample 3000, drop contaminants) ----------
ec <- readRDS("D:/escc_sc/out/endo_annotated.rds")
ec <- subset(ec, subset = !grepl("^Contam", as.character(celltype)))
ec$grp <- as.character(ec$celltype)
set.seed(1)
idx <- unlist(lapply(split(colnames(ec), ec$grp), function(x) sample(x, min(length(x), 350))))
ec <- subset(ec, cells = idx)
cat("EC cells:", ncol(ec), "\n"); print(table(ec$grp))

# ---------- T cells (8000 random cells, full genes) ----------
tc <- readMM("D:/escc_sc/cellchat/tcell_subset.mtx")
tg <- readLines("D:/escc_sc/cellchat/tcell_genes.txt")
tb <- readLines("D:/escc_sc/cellchat/tcell_barcodes.txt")
rownames(tc) <- tg; colnames(tc) <- tb
tc <- as(tc, "CsparseMatrix")

# T subset labels. Default is now an explicit unclassified category.
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
lab[cd8 & tex] <- "CD8 exhausted"   # >=2 exhaustion-associated markers
lab[treg] <- "Treg"
lab[nk] <- "NK-like"
names(lab) <- colnames(tc)
print(table(lab))
write.csv(as.data.frame(table(lab)), file.path(OUT, "tcell_group_counts_v3_5.csv"),
          row.names = FALSE)

# ---------- merge ----------
ec_counts <- LayerData(ec, assay = "RNA", layer = "counts")
common <- intersect(rownames(ec_counts), rownames(tc))
cat("common genes:", length(common), "\n")
mat <- cbind(ec_counts[common, ], tc[common, ])
grp <- c(ec$grp, lab[colnames(tc)])
meta <- data.frame(labels = grp, samples = c(as.character(ec$sample), sub("-.*$", "", colnames(tc))),
                   row.names = colnames(mat))
cat("total cells:", ncol(mat), "\n")
print(table(meta$labels))

obj <- CreateSeuratObject(counts = mat, meta.data = meta)
obj <- NormalizeData(obj, verbose = FALSE)
data.input <- LayerData(obj, assay = "RNA", layer = "data")
rm(mat, obj); gc()

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 50, seed.use = 1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(OUT, "cellchat_net.csv"), row.names = FALSE)
cat("retained ligand-receptor pairs:", nrow(df.net), "\n")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(OUT, "cellchat.rds"))

# ---------- figures (intermediate, not submitted) ----------
png(file.path(OUT, "net_circle_count.png"), width = 2000, height = 1800, res = 300)
par(mfrow = c(1, 1), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of retained ligand-receptor pairs (EC -> T)")
dev.off()

png(file.path(OUT, "net_circle_weight.png"), width = 2000, height = 1800, res = 300)
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Aggregate CellChat score (EC -> T, arbitrary units)")
dev.off()

p <- netVisual_heatmap(cellchat, measure = "weight")
ggsave(file.path(OUT, "heatmap_weight.png"), p, width = 7, height = 5.5, dpi = 300)

p2 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(file.path(OUT, "role_scatter.png"), p2, width = 6.5, height = 5, dpi = 300)

rk <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = FALSE,
              return.data = TRUE)
write.csv(as.data.frame(rk$signaling.contribution),
          file.path(OUT, "pathway_contribution.csv"), row.names = FALSE)
p3 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = FALSE)
ggsave(file.path(OUT, "pathway_rank.png"), p3, width = 7, height = 4.5, dpi = 300)

cat("\nDONE CellChat v3.5\n")
