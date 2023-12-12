## Description: Draw plots showing DvP scores on PCA plot and boxplot with conditions

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)

figureDir <- file.path("figures", "manuscript")
vsd <- readRDS("data/processed/deseq/limma_batch_normalized.rds")
pca_df <- pca_df <- plotPCA(vsd, returnData=T)
esdf <- read_csv("data/gsea-results/Hallmark_Reactome_DvP_GSVA_Scores.csv")
esdf <- esdf %>% 
  inner_join(pca_df, by = c("sampleID" = "name", "condition")) %>%
  # filter(condition != "AK_IEC", condition != "AK_IEC_SCC") %>%
  mutate(condition = factor(condition, levels = c("NS", "AK", "AK_IEC", "AK_IEC_SCC", "IEC", "SCC", "KA")))

condition_colors <- readRDS("data/condition_color_palette.rds")

pca_plot <- esdf %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = condition), size = 3, pch=21) +
  theme_classic() +
  scale_fill_manual(values = condition_colors) +
  labs(fill = "")

pca_score_plot <- esdf %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = DvP_Score), size = 3, pch=21) +
  theme_classic() +
  scale_fill_viridis() +
  labs(fill = "DvP Score") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  )


dvpBoxplot <- esdf %>%
  ggplot(aes(condition, DvP_Score)) +
  geom_boxplot(aes(fill = condition)) +
  theme_classic() +
  labs(x = "", y = "DvP Score") +
  guides(fill = "none") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = condition_colors)
dvpBoxplot
lm(DvP_Score ~ condition, data = esdf) %>% summary()
kruskal.test(DvP_Score ~ condition, data = esdf %>% filter(condition != "AK_IEC", condition != "AK_IEC_SCC"))

# esdf %>%
#   ggplot(aes(condition, DvP_Score)) +
#   geom_violin(aes(fill = condition)) +
#   geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
#   theme_classic() +
#   labs(x = "", y = "DvP Score", title = "DvP Score") +
#   guides(fill = "none") +
#   theme(plot.title = element_text(hjust = .5)) +
#   scale_fill_manual(values = condition_colors)
# lm(DvP_Score ~ condition, data = esdf) %>% summary()
# kruskal.test(DvP_Score ~ condition, data = esdf %>% filter(condition != "AK_IEC", condition != "AK_IEC_SCC"))

ggsave(
  filename = file.path(figureDir, "Figure1_DvP_PCA_Plot.svg"),
  plot = pca_score_plot,
  width = 5,
  height = 3
)

