# Pseudobulk DE sensitivity: log2-CPM limma (primary) vs TMM + voom precision weights.
# All English; run with Rscript.exe. Inputs are the archived analysis objects.
suppressPackageStartupMessages({
  library(Seurat); library(limma); library(edgeR); library(dplyr); library(tidyr)
  library(clusterProfiler); library(org.Hs.eg.db)
})

OUT <- "D:/escc_sc/out_v43"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logp <- function(x) cat(format(Sys.time(), "%H:%M:%S"), x, "\n",
                        file = file.path(OUT, "progress.log"), append = TRUE)
logp("START")
obj <- readRDS("D:/escc_sc/out/endo_annotated.rds")
obj <- subset(obj, subset = type == "Tumor")
obj$sc <- paste(obj$sample, as.character(obj$celltype), sep = "|")
tab <- as.data.frame(table(obj$sc)); colnames(tab) <- c("sc", "n")
keep_sc <- tab$sc[tab$n >= 10]
pb <- AggregateExpression(obj, group.by = "sc", assays = "RNA", slot = "counts")$RNA
colnames(pb) <- gsub("-", "_", colnames(pb))
pb <- pb[, colnames(pb) %in% keep_sc, drop = FALSE]
meta <- data.frame(sc = colnames(pb), stringsAsFactors = FALSE) %>%
  separate(sc, into = c("sample", "celltype"), sep = "\\|", extra = "merge") %>%
  mutate(celltype = gsub("-", "_", celltype))
rownames(meta) <- meta$sc
cat("pseudobulk matrix:", dim(pb), "\n")
logp("pseudobulk built")

ref487 <- read.csv("D:/escc_sc/out/Tip_signature_shared.csv", stringsAsFactors = FALSE)$gene

run_one <- function(a, b, tag) {
  pat <- meta %>% filter(celltype %in% c(a, b)) %>%
    group_by(sample) %>% filter(n() == 2) %>% pull(sample) %>% unique()
  if (length(pat) < 8) return(NULL)
  sel <- meta$sample %in% pat & meta$celltype %in% c(a, b)
  m <- pb[, sel, drop = FALSE]
  md <- meta[sel, ]
  md$group <- factor(md$celltype, levels = c(b, a))
  design <- model.matrix(~ patient + group, data = transform(md, patient = factor(md$sample)))
  group_col <- grep("^group", colnames(design), value = TRUE)[1]
  cpm_all <- t(t(m) / colSums(m)) * 1e6
  keep <- rowSums(cpm_all > 1) >= 0.25 * ncol(cpm_all)
  m <- as.matrix(m[keep, , drop = FALSE])
  cat(tag, "pairs =", length(pat), "genes =", nrow(m), "\n")
  logp(paste(tag, "start, genes =", nrow(m)))

  # Method 1: primary (log2-CPM + ordinary limma + eBayes)
  lcpm <- log2(as.matrix(cpm_all[keep, , drop = FALSE]) + 1)
  f1 <- eBayes(lmFit(lcpm, design))
  t1 <- topTable(f1, coef = group_col, number = Inf, sort.by = "P")
  t1$gene <- rownames(t1)

  logp(paste(tag, "baseline done"))
  # Method 2: TMM + voom + limma
  dge <- calcNormFactors(DGEList(counts = m), method = "TMM")
  logp(paste(tag, "TMM done"))
  v <- voom(dge, design)
  logp(paste(tag, "voom transform done"))
  f2 <- eBayes(lmFit(v, design))
  t2 <- topTable(f2, coef = group_col, number = Inf, sort.by = "P")
  t2$gene <- rownames(t2)

  logp(paste(tag, "voom done"))

  get_up <- function(tt, pcol, fccol) {
    ok <- !is.na(tt[[pcol]]) & !is.na(tt[[fccol]]) &
      tt[[pcol]] < 0.05 & tt[[fccol]] > 1
    ok
  }
  list(
    tag = tag, n_genes = nrow(m), n_pairs = length(pat),
    baseline = list(tab = t1, up = t1$gene[get_up(t1, "adj.P.Val", "logFC")]),
    voom = list(tab = t2, up = t2$gene[get_up(t2, "adj.P.Val", "logFC")])
  )
}

contrasts <- list(
  c("Tip_angiogenic", "Venous", "Tip_vs_Venous"),
  c("Tip_angiogenic", "Capillary", "Tip_vs_Capillary")
)
res <- list()
for (cc in contrasts) {
  r <- run_one(cc[1], cc[2], cc[3])
  if (!is.null(r)) res[[cc[3]]] <- r
}

methods <- c("baseline", "voom")
summary_rows <- tip_rows <- list()
sig_by_method <- setNames(vector("list", length(methods)), methods)
for (me in methods) {
  u1 <- res$Tip_vs_Venous[[me]]$up
  u2 <- res$Tip_vs_Capillary[[me]]$up
  tip <- intersect(u1, u2)
  sig_by_method[[me]] <- tip
  ov <- length(intersect(tip, ref487))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    method = me, contrast = "TIP_SIGNATURE", n_genes_tested = NA, n_pairs = NA,
    n_sig_up = NA, n_tip_signature = length(tip), overlap_with_primary_487 = ov,
    pct_of_primary_retained = round(100 * ov / length(ref487), 1),
    jaccard_vs_primary = round(length(intersect(tip, ref487)) /
                                 length(union(tip, ref487)), 4))
  for (cc in contrasts) {
    r <- res[[cc[3]]]
    if (is.null(r)) next
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      method = me, contrast = cc[3], n_genes_tested = r$n_genes,
      n_pairs = r$n_pairs, n_sig_up = length(r[[me]]$up),
      n_tip_signature = NA, overlap_with_primary_487 = NA,
      pct_of_primary_retained = NA, jaccard_vs_primary = NA)
  }
  for (g in tip) tip_rows[[length(tip_rows) + 1]] <- data.frame(method = me, gene = g)
}

# GO BP enrichment for each tip signature (same settings as primary)
go_rows <- list()
for (me in methods) {
  g <- sig_by_method[[me]]
  eg <- bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 readable = TRUE)
  gd <- as.data.frame(go)
  if (nrow(gd) == 0) next
  for (i in seq_len(min(8, nrow(gd)))) {
    go_rows[[length(go_rows) + 1]] <- data.frame(
      method = me, rank = i, Description = gd$Description[i],
      GeneRatio = gd$GeneRatio[i], p.adjust = signif(gd$p.adjust[i], 3),
      Count = gd$Count[i])
  }
}

write.csv(do.call(rbind, summary_rows), file.path(OUT, "voom_sensitivity_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, tip_rows), file.path(OUT, "voom_sensitivity_tip_genes.csv"), row.names = FALSE)
write.csv(do.call(rbind, go_rows), file.path(OUT, "voom_sensitivity_GO.csv"), row.names = FALSE)
cat("\n=== SUMMARY ===\n")
print(do.call(rbind, summary_rows), row.names = FALSE)
cat("\n=== GO top terms ===\n")
print(do.call(rbind, go_rows), row.names = FALSE)
cat("\nDONE\n")
