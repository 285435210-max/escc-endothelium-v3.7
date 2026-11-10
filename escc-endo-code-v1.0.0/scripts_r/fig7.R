# Figure 5 (script Fig7_communication): endothelial-T cell communication - v3.5 layout
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(stringr)
})
FIG <- "D:/escc_sc/figs"; CC <- "D:/escc_sc/cellchat_out"
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0))

net <- read.csv(file.path(CC, "cellchat_net.csv"), check.names = FALSE)
ec <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
        "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic")
tt <- c("CD8 T","CD8 exhausted","CD4 T","Treg","NK-like")
DISP <- c(Arterial = "Arterial", Capillary = "Capillary", Venous = "Venous",
          Venous_inflamed = "Inflamed venule", Tip_angiogenic = "Tip-like",
          Matrix_EC = "Matrix-remodeling", IFNgamma_EC = "Interferon-responsive",
          Proliferating = "Proliferating", Lymphatic = "Lymphatic")
sub <- net %>% filter(source %in% ec, target %in% tt)

# ---- a: heatmap (count of retained ligand-receptor pairs) ----
agg <- sub %>% group_by(source, target) %>%
  summarise(n = n(), prob = sum(prob), .groups = "drop") %>%
  mutate(source = factor(unname(DISP[source]), levels = rev(unname(DISP[ec]))),
         target = factor(target, levels = tt))
p7a <- ggplot(agg, aes(target, source, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n), size = 1.9, colour = "grey10") +
  scale_fill_gradient(low = "#fdf0ea", high = "#8c2d19") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = "Retained pairs",
       title = "Number of retained ligand-receptor pairs") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(),
             axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
             axis.text.y = element_text(size = 6.2),
             legend.key.height = unit(4, "mm"), legend.key.width = unit(2.2, "mm"),
             legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
             aspect.ratio = 0.85)

# ---- b: pathway ranking ----
pw <- sub %>% group_by(pathway_name) %>%
  summarise(prob = sum(prob), n = n(), .groups = "drop") %>% arrange(desc(prob)) %>% head(11)
pw$pathway_name <- factor(pw$pathway_name, levels = rev(pw$pathway_name))
p7b <- ggplot(pw, aes(prob, pathway_name)) +
  geom_col(fill = "#4a6ca8", colour = "black", linewidth = 0.15, width = 0.7) +
  geom_text(aes(label = paste0("n = ", n)), hjust = -0.25, size = 1.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(x = "Aggregate CellChat score (arbitrary units)", y = NULL,
       title = "Signaling pathways") +
  th + theme(axis.text.y = element_text(size = 6.2))

# ---- c: top LR pairs to CD8 T cells, coloured by pathway ----
cd8 <- sub %>% filter(target %in% c("CD8 T", "CD8 exhausted")) %>%
  group_by(interaction_name_2, pathway_name) %>%
  summarise(prob = sum(prob), .groups = "drop") %>%
  arrange(desc(prob)) %>% head(12)
cd8$lab <- gsub("_", " - ", cd8$interaction_name_2)
cd8$lab <- factor(cd8$lab, levels = rev(cd8$lab))
pw_cols <- c("MHC-I" = "#b03a2e", "MIF" = "#7b4fa0", "APP" = "#d98b5f",
             "COLLAGEN" = "#4a9bc4", "GALECTIN" = "#2e7d32", "CD99" = "#8c6d31",
             "MHC-II" = "#c0553f")
p7c <- ggplot(cd8, aes(prob, lab, fill = pathway_name)) +
  geom_col(colour = "black", linewidth = 0.15, width = 0.72) +
  scale_fill_manual(values = pw_cols, na.value = "grey70") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Aggregate CellChat score (arbitrary units)", y = NULL, fill = "Pathway",
       title = "Ligand-receptor pairs targeting CD8 T cells") +
  th + theme(axis.text.y = element_text(size = 5.8),
             legend.position = "right", legend.key.size = unit(2.5, "mm"),
             legend.title = element_text(size = 6), legend.text = element_text(size = 5.5))

fig7 <- (p7a | p7b) / p7c +
  plot_layout(widths = c(1.15, 1), heights = c(1, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Fig7_communication.pdf"), fig7, width = 180, height = 150, units = "mm")
ggsave(file.path(FIG, "Fig7_communication.png"), fig7, width = 180, height = 150, units = "mm", dpi = 600)
cat("Fig7 done\n")
