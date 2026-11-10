suppressPackageStartupMessages({library(CellChat); library(ggplot2); library(ComplexHeatmap); library(dplyr)})
OUT <- "D:/escc_sc/cellchat_out"
cc <- readRDS(file.path(OUT, "cellchat.rds"))
cat("groups:", length(levels(cc@idents)), " interactions:", nrow(subsetCommunication(cc)), "\n")

png(file.path(OUT, "heatmap_weight.png"), width = 2400, height = 2000, res = 300)
print(draw(netVisual_heatmap(cc, measure = "weight", color.heatmap = "Reds")))
dev.off()

png(file.path(OUT, "heatmap_count.png"), width = 2400, height = 2000, res = 300)
print(draw(netVisual_heatmap(cc, measure = "count", color.heatmap = "Reds")))
dev.off()

cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
saveRDS(cc, file.path(OUT, "cellchat.rds"))
p <- netAnalysis_signalingRole_scatter(cc, title = "Signaling roles: endothelium vs T cells")
ggsave(file.path(OUT, "role_scatter.png"), p, width = 6.5, height = 5.2, dpi = 300)

# EC -> T 的 pathway 汇总
net <- subsetCommunication(cc)
ec_groups <- c("Arterial","Capillary","Venous","Venous_inflamed","Tip_angiogenic",
               "Matrix_EC","IFNgamma_EC","Proliferating","Lymphatic")
t_groups  <- c("CD8 T","CD8 exhausted","CD4 T","Treg","NK-like")
sub <- net %>% filter(source %in% ec_groups, target %in% t_groups)
agg <- sub %>% group_by(pathway_name) %>%
  summarise(n_pairs = n(), total_prob = sum(prob), total_pval = mean(pval), .groups="drop") %>%
  arrange(desc(total_prob))
write.csv(agg, file.path(OUT, "pathway_EC_to_T.csv"), row.names = FALSE)
cat("\nTop pathways (EC -> T):\n"); print(head(agg, 10))
p4 <- ggplot(head(agg, 15), aes(x = reorder(pathway_name, total_prob), y = total_prob)) +
  geom_col(fill = "#4a6fa5") + coord_flip() + theme_bw() +
  labs(x = NULL, y = "Total communication probability", title = "EC -> T signaling pathways")
ggsave(file.path(OUT, "pathway_bar.png"), p4, width = 6.8, height = 4.6, dpi = 300)

lr <- sub %>% group_by(interaction_name) %>%
  summarise(n = n(), prob = sum(prob), .groups="drop") %>% arrange(desc(prob))
write.csv(head(lr, 30), file.path(OUT, "top_LR_EC_to_T.csv"), row.names = FALSE)
cat("\nTop LR pairs:/n"); print(head(lr, 12))

# 关键通路各亚群贡献热图
key <- head(agg$pathway_name, 8)
png(file.path(OUT, "heatmap_key_pathways.png"), width = 2600, height = 1600, res = 300)
print(draw(netVisual_heatmap(cc, measure = "weight", signaling = key[1:4], color.heatmap = "Reds")))
dev.off()
cat("\nDONE figs\n")
