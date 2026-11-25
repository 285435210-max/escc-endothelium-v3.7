# GO enrichment sensitivity: genome-wide default background vs tested-gene universe.
# All ASCII; run with Rscript.exe. Derives universe from the archived pseudobulk objects.
suppressPackageStartupMessages({
  library(Seurat); library(limma); library(edgeR); library(dplyr); library(tidyr)
  library(clusterProfiler); library(org.Hs.eg.db)
})
OUT <- "D:/escc_sc/out_v43"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(OUT, "go_universe_progress.log")
logp <- function(x) cat(format(Sys.time(), "%H:%M:%S"), x, "\n", file = logf, append = TRUE)
logp("START go tested-universe")

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
logp(paste("pseudobulk built", nrow(pb), ncol(pb)))

tested_lists <- list()
for (cc in list(c("Tip_angiogenic", "Venous", "Tip_vs_Venous"),
                c("Tip_angiogenic", "Capillary", "Tip_vs_Capillary"))) {
  a <- cc[1]; b <- cc[2]; tag <- cc[3]
  pat <- meta %>% filter(celltype %in% c(a, b)) %>%
    group_by(sample) %>% filter(n() == 2) %>% pull(sample) %>% unique()
  sel <- meta$sample %in% pat & meta$celltype %in% c(a, b)
  m <- pb[, sel, drop = FALSE]; md <- meta[sel, ]
  md$group <- factor(md$celltype, levels = c(b, a))
  cpm_all <- t(t(m) / colSums(m)) * 1e6
  keep <- rowSums(cpm_all > 1) >= 0.25 * ncol(cpm_all)
  genes <- rownames(m)[keep]
  tested_lists[[tag]] <- genes
  logp(paste(tag, "pairs", length(pat), "tested genes", length(genes)))
}
univ <- Reduce(intersect, tested_lists)
logp(paste("common tested universe", length(univ)))
write.csv(data.frame(gene = univ), file.path(OUT, "go_tested_universe_genes.csv"), row.names = FALSE)

ref487 <- read.csv("D:/escc_sc/out/Tip_signature_shared.csv", stringsAsFactors = FALSE)$gene
voom241 <- read.csv(file.path(OUT, "voom_sensitivity_tip_genes.csv"), stringsAsFactors = FALSE)
voom241 <- unique(voom241$gene[voom241$method == "voom"])

run_go <- function(genes, label) {
  gl <- intersect(genes, univ)
  eg <- bitr(gl, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  eg_univ <- bitr(univ, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  r <- enrichGO(gene = eg$ENTREZID, universe = eg_univ$ENTREZID,
                OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  d <- as.data.frame(r)
  if (nrow(d) == 0) d <- data.frame(ID = NA, Description = NA, p.adjust = NA)
  d$gene_list <- label
  d$n_input <- length(gl)
  d$n_universe <- length(univ)
  return(d)
}

d1 <- run_go(ref487, "baseline_487")
logp(paste("baseline GO terms", nrow(d1)))
d2 <- run_go(voom241, "voom_241")
logp(paste("voom GO terms", nrow(d2)))
res <- rbind(d1, d2)
write.csv(res, file.path(OUT, "go_tested_universe_result.csv"), row.names = FALSE)
cat("=== baseline 487 top terms ===\n")
print(head(d1[, c("Description", "p.adjust", "Count")], 6))
cat("=== voom 241 top terms ===\n")
print(head(d2[, c("Description", "p.adjust", "Count")], 6))
logp("DONE")
