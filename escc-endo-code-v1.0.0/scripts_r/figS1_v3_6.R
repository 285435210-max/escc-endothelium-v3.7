# Supplementary Figure S1, v3.6: display names instead of internal code names.
# Replots from cached SingleR / CellMarker result tables (no rerun of SingleR).
suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(patchwork); library(dplyr)
})
OUT <- "D:/escc_sc/annot_check"
pct <- as.matrix(read.csv(file.path(OUT, "singler_rowpercent.csv"),
                          row.names = 1, check.names = FALSE))
z <- as.matrix(read.csv(file.path(OUT, "cellmarker_zscores.csv"),
                        row.names = 1, check.names = FALSE))

disp <- c(Arterial = "Arterial", Capillary = "Capillary", Venous = "Venous",
          Venous_inflamed = "Inflamed venule", Tip_angiogenic = "Tip-like",
          Matrix_EC = "Matrix-remodeling", IFNgamma_EC = "Interferon-responsive",
          Proliferating = "Proliferating", Lymphatic = "Lymphatic",
          Contam_pericyte = "Excluded: pericyte",
          Contam_fibroblast = "Excluded: fibroblast",
          Contam_immune = "Excluded: immune")
ord <- c("Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
         "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic",
         "Contam_pericyte", "Contam_fibroblast", "Contam_immune")
ord <- ord[ord %in% rownames(pct)]
ord_disp <- unname(disp[ord])
P <- pct[ord, , drop = FALSE]
rownames(P) <- ord_disp
z2 <- z[ord, , drop = FALSE]
rownames(z2) <- ord_disp

th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 8, face = "bold"),
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5))

keep <- colnames(P)[apply(P, 2, max) >= 2]
rest <- setdiff(colnames(P), keep)
dfA <- as.data.frame(P[, keep, drop = FALSE]); dfA$state <- rownames(dfA)
if (length(rest)) dfA$Other <- rowSums(P[, rest, drop = FALSE])
long <- pivot_longer(dfA, -state, names_to = "ref", values_to = "pct")
long$state <- factor(long$state, levels = rev(ord_disp))
long$ref <- factor(long$ref, levels = c(keep, if (length(rest)) "Other"))
pal <- c(`Endothelial cells` = "#1b6ca8", Adipocytes = "#e3c08d", Fibroblasts = "#a9775a",
         HSC = "#7b4fa0", `CD8+ T-cells` = "#4a9b6e", Other = "#c9c9c9")

pA <- ggplot(long, aes(pct, state, fill = ref)) +
  geom_col(width = 0.72) +
  geom_text(data = data.frame(state = factor(ord_disp, levels = rev(ord_disp)),
                              lab = sprintf("%.0f", P[, "Endothelial cells"])),
            aes(x = 104, y = state, label = lab), inherit.aes = FALSE, size = 2.1, hjust = 0) +
  geom_hline(yintercept = 3.5, linewidth = 0.3, colour = "grey40") +
  scale_fill_manual(values = pal, name = NULL,
                    breaks = c("Endothelial cells", "Adipocytes", "Fibroblasts", "HSC",
                               "CD8+ T-cells", "Other")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 116), breaks = c(0, 50, 100)) +
  labs(x = "Cells assigned to each reference label (%)", y = NULL,
       title = "SingleR vs Blueprint ENCODE") +
  th + theme(legend.position = "right", legend.key.height = unit(2.6, "mm"))

dfB <- as.data.frame(z2); dfB$state <- rownames(dfB)
longB <- pivot_longer(dfB, -state, names_to = "set", values_to = "z")
longB$state <- factor(longB$state, levels = rev(ord_disp))
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
ggsave(file.path(OUT, "FigS1_annotation_validation.png"), fig,
       width = 180, height = 120, units = "mm", dpi = 600)
ggsave(file.path(OUT, "FigS1_annotation_validation.pdf"), fig,
       width = 180, height = 120, units = "mm")
cat("DONE FigS1 v3.6\n")
