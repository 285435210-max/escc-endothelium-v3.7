# Pseudobulk differential expression between endothelial states (paired by patient)
suppressPackageStartupMessages({
  library(Seurat); library(limma); library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "D:/escc_sc/out"
obj <- readRDS(file.path(OUT, "endo_annotated.rds"))
obj <- subset(obj, subset = type == "Tumor")
obj$sc <- paste(obj$sample, as.character(obj$celltype), sep = "|")

# cells per sample x celltype, keep >= 10
tab <- as.data.frame(table(obj$sc))
colnames(tab) <- c("sc", "n")
keep_sc <- tab$sc[tab$n >= 10]
cat("pseudobulk groups (>=10 cells):", length(keep_sc), "\n")

pb <- AggregateExpression(obj, group.by = "sc", assays = "RNA", slot = "counts")$RNA
colnames(pb) <- gsub("-", "_", colnames(pb))          # Seurat sanitizes underscores in identities
pb <- pb[, colnames(pb) %in% keep_sc, drop = FALSE]
meta <- data.frame(sc = colnames(pb), stringsAsFactors = FALSE) %>%
  separate(sc, into = c("sample", "celltype"), sep = "\\|", extra = "merge") %>%
  mutate(celltype = gsub("-", "_", celltype))
rownames(meta) <- meta$sc
cat("pseudobulk matrix:", dim(pb), "\n")
write.csv(as.data.frame(table(meta$celltype)), file.path(OUT, "pseudobulk_group_sizes.csv"), row.names = FALSE)

run_contrast <- function(a, b, tag) {
  pat <- meta %>% filter(celltype %in% c(a, b)) %>%
    group_by(sample) %>% filter(n() == 2) %>% pull(sample) %>% unique()
  if (length(pat) < 8) { cat("skip", tag, ": paired patients =", length(pat), "\n"); return(NULL) }
  sel <- meta$sample %in% pat & meta$celltype %in% c(a, b)
  m <- pb[, sel, drop = FALSE]
  md <- meta[sel, ]
  md$group <- factor(md$celltype, levels = c(b, a))
  cat("\n", tag, ": pairs =", length(pat), " groups =", paste(table(md$group), collapse = "/"), "\n")

  cpm <- t(t(m) / colSums(m)) * 1e6
  lcpm <- log2(cpm + 1)
  keep <- rowSums(cpm > 1) >= 0.25 * ncol(cpm)
  lcpm <- lcpm[keep, , drop = FALSE]
  cat("  genes tested:", nrow(lcpm), "\n")

  design <- model.matrix(~ patient + group, data = transform(md, patient = factor(md$sample)))
  fit <- lmFit(lcpm, design)
  fit <- eBayes(fit)
  coef_name <- paste0("group", a)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  write.csv(tt, file.path(OUT, paste0("DE_", tag, ".csv")), row.names = FALSE)
  up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 1)
  dn <- sum(tt$adj.P.Val < 0.05 & tt$logFC < -1)
  cat("  DE (FDR<0.05, |logFC|>1): up =", up, " down =", dn, "\n")

  p <- tt %>% mutate(sig = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, ifelse(logFC > 0, "Up", "Down"), "NS"))
  g <- ggplot(p, aes(logFC, -log10(P.Value), color = sig)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_color_manual(values = c(Up = "#b03a2e", Down = "#2e6da4", NS = "grey80")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", size = 0.3) +
    labs(title = paste0(tag, " (n=", length(pat), " paired patients)"),
         x = "log2 fold change", y = "-log10 P") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(OUT, paste0("Volcano_", tag, ".png")), g, width = 5, height = 4.2, dpi = 300)
  tt
}

res <- list()
res[["Tip_vs_Venous"]]      <- run_contrast("Tip_angiogenic", "Venous", "Tip_vs_Venous")
res[["Matrix_vs_Venous"]]   <- run_contrast("Matrix_EC", "Venous", "Matrix_vs_Venous")
res[["IFNg_vs_Venous"]]     <- run_contrast("IFNgamma_EC", "Venous", "IFNg_vs_Venous")
res[["Tip_vs_Capillary"]]   <- run_contrast("Tip_angiogenic", "Capillary", "Tip_vs_Capillary")

# shared tip signature (robust across two references)
if (!is.null(res$Tip_vs_Venous) && !is.null(res$Tip_vs_Capillary)) {
  a <- res$Tip_vs_Venous; b <- res$Tip_vs_Capillary
  shared <- merge(a[, c("gene", "logFC", "adj.P.Val")], b[, c("gene", "logFC", "adj.P.Val")],
                  by = "gene", suffixes = c("_vsVenous", "_vsCap"))
  sig <- shared %>% filter(adj.P.Val_vsVenous < 0.05, adj.P.Val_vsCap < 0.05,
                           logFC_vsVenous > 1, logFC_vsCap > 1) %>% arrange(desc(logFC_vsVenous))
  write.csv(sig, file.path(OUT, "Tip_signature_shared.csv"), row.names = FALSE)
  cat("\nshared tip-up genes (FDR<0.05 both contrasts):", nrow(sig), "\n")
  print(head(sig, 20))
}

cat("\nDONE DE\n")
