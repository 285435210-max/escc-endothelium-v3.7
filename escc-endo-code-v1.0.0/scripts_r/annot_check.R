# Reference-based annotation check: SingleR (Blueprint ENCODE) + CellMarker projection
suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(ComplexHeatmap); library(circlize); library(dplyr)
})
set.seed(1)
OUT <- "D:/escc_sc/annot_check"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS("D:/escc_sc/out/endo_annotated.rds")
cat("cells:", ncol(obj), " labels:", length(unique(obj$celltype)), "\n")

mat <- GetAssayData(obj, assay = "RNA", layer = "data")

# ---------- panel A: SingleR vs Blueprint ENCODE ----------
ref <- BlueprintEncodeData()
cat("reference:", dim(ref), "\n")
pred <- SingleR(test = mat, ref = ref, labels = ref$label.main, assay.type.test = 1)
obj$singleR <- pred$labels

tab <- table(our = obj$celltype, ref = pred$labels)
pct <- round(100 * prop.table(tab, 1), 1)
write.csv(as.data.frame.matrix(tab), file.path(OUT, "singler_counts.csv"))
write.csv(as.data.frame.matrix(pct), file.path(OUT, "singler_rowpercent.csv"))

ec_tab <- data.frame(state = rownames(pct),
                     pct_endothelial = pct[, "Endothelial cells"],
                     n = as.integer(rowSums(tab)))
write.csv(ec_tab, file.path(OUT, "EC_call_rate.csv"), row.names = FALSE)
cat("\n--- endothelial call rate (Blueprint) ---\n"); print(ec_tab)
w <- sum(tab[, "Endothelial cells"]) / sum(tab)
cat(sprintf("weighted overall: %.1f%% of %d cells\n", 100 * w, sum(tab)))

top <- t(sapply(rownames(pct), function(r) {
  i <- which.max(pct[r, ])
  c(label = colnames(pct)[i], pct = pct[r, i])
}))
cat("\n--- top SingleR label per state ---\n"); print(top)

# ---------- panel B: CellMarker projection ----------
read_gmt <- function(f) {
  ln <- strsplit(readLines(f, warn = FALSE), "\t")
  setNames(lapply(ln, function(x) unique(x[-c(1, 2)])), sapply(ln, `[`, 1))
}
cm <- read_gmt("D:/escc_sc/data/CellMarker_2024.gmt")
pick <- c(Arterial = "Endothelial Cell Artery Human",
          Venous = "Endothelial Cell Vein Human",
          Capillary = "Endothelial Cell Capillary Human",
          Lymphatic = "Endothelial Cell Lymphatic Vessel Human",
          Esophagus_EC = "Endothelial Cell Esophagus Human",
          Fibroblast = "Fibroblast Esophagus Human",
          Pericyte = "Pericyte Stomach Human",
          `T cell` = "CD8+ T Cell Blood Human",
          Macrophage = "Macrophage Blood Human")
sets <- lapply(pick, function(k) intersect(cm[[k]], rownames(obj)))
cat("\nset sizes:", sapply(sets, length), "\n")

obj <- AddModuleScore(obj, features = sets, name = "ref_", seed = 1)
sc <- obj@meta.data[, paste0("ref_", seq_along(sets)), drop = FALSE]
colnames(sc) <- names(sets)
sc$celltype <- obj$celltype
avg <- sc %>% group_by(celltype) %>% summarise(across(everything(), mean), .groups = "drop")
m <- as.data.frame(avg[, -1]); rownames(m) <- avg$celltype
z <- scale(m)                      # z-score per reference set
write.csv(round(as.matrix(m), 3), file.path(OUT, "cellmarker_scores.csv"))
write.csv(round(as.matrix(z), 2), file.path(OUT, "cellmarker_zscores.csv"))

# ---------- figure: ggplot version ----------
suppressPackageStartupMessages({library(ggplot2); library(tidyr); library(patchwork)})
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 8, face = "bold"),
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5))

ord <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic","Matrix_EC",
         "IFNgamma_EC","Proliferating","Lymphatic",
         "Contam_pericyte","Contam_fibroblast","Contam_immune")
ord <- ord[ord %in% rownames(pct)]
P <- matrix(as.numeric(pct), nrow = nrow(pct), dimnames = dimnames(pct))
P <- P[ord, , drop = FALSE]

keep <- colnames(P)[apply(P, 2, max) >= 2]
rest <- setdiff(colnames(P), keep)
dfA <- as.data.frame(P[, keep, drop = FALSE]); dfA$state <- rownames(dfA)
if (length(rest)) dfA$Other <- rowSums(P[, rest, drop = FALSE])
long <- pivot_longer(dfA, -state, names_to = "ref", values_to = "pct")
long$state <- factor(long$state, levels = rev(ord))
long$ref <- factor(long$ref, levels = c(keep, if (length(rest)) "Other"))
pal <- c(`Endothelial cells` = "#1b6ca8", Adipocytes = "#e3c08d", Fibroblasts = "#a9775a",
         HSC = "#7b4fa0", `CD8+ T-cells` = "#4a9b6e", Other = "#c9c9c9")

pA <- ggplot(long, aes(pct, state, fill = ref)) +
  geom_col(width = 0.72) +
  geom_text(data = data.frame(state = factor(ord, levels = rev(ord)),
                              lab = sprintf("%.0f", pct[ord, "Endothelial cells"])),
            aes(x = 104, y = state, label = lab), inherit.aes = FALSE, size = 2.1, hjust = 0) +
  geom_hline(yintercept = 3.5, linewidth = 0.3, colour = "grey40") +
  scale_fill_manual(values = pal, name = NULL,
                    breaks = c("Endothelial cells", "Adipocytes", "Fibroblasts", "HSC", "CD8+ T-cells", "Other")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 116), breaks = c(0, 50, 100)) +
  labs(x = "Cells assigned to each reference label (%)", y = NULL,
       title = "SingleR vs Blueprint ENCODE") +
  th + theme(legend.position = "right", legend.key.height = unit(2.6, "mm"))

z2 <- z[ord, , drop = FALSE]
dfB <- as.data.frame(z2); dfB$state <- rownames(dfB)
longB <- pivot_longer(dfB, -state, names_to = "set", values_to = "z")
longB$state <- factor(longB$state, levels = rev(ord))
longB$set <- factor(longB$set, levels = colnames(z2))

pB <- ggplot(longB, aes(set, state, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(abs(z) >= 1.5, sprintf("%.1f", z), ""),
                colour = abs(z) >= 1.9), size = 1.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  geom_hline(yintercept = 3.5, linewidth = 0.3, colour = "grey40") +
  scale_fill_gradient2(low = "#2c5f8a", mid = "white", high = "#b03a2e", midpoint = 0,
                       limits = c(-3, 3), name = "z-score") +
  labs(x = NULL, y = NULL, title = "CellMarker projection") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
             axis.line = element_blank(), axis.ticks = element_blank())

fig <- pA / pB + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "FigS1_annotation_validation.png"), fig, width = 180, height = 120, units = "mm", dpi = 600)
ggsave(file.path(OUT, "FigS1_annotation_validation.pdf"), fig, width = 180, height = 120, units = "mm")
cat("\nDONE\n")
