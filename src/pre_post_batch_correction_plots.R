## Description: Draw PCA plots pre and post limma batch correction


library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(stringr)

theme_set(theme_classic())

dataDir <- file.path("data", "processed", "deseq")
figureDir <- file.path("figures", "manuscript")
pre_vsd <- readRDS(file.path(dataDir, "vst_normalized_counts.rds"))
post_vsd <- readRDS(file.path(dataDir, "limma_batch_normalized.rds"))

condition_colors <- readRDS("data/condition_color_palette.rds")

pre_df <- plotPCA(pre_vsd, intgroup = c("study", "condition"), returnData = TRUE)
post_df <- plotPCA(post_vsd, intgroup = c("study", "condition"), returnData = TRUE)

pre_df$study <- str_replace(pre_df$study, "_", " ")
post_df$study <- str_replace(post_df$study, "_", " ")

p1 <- pre_df %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = study), size=3, pch=21) +
  labs(fill = "", title = "Study") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  )
p2 <- pre_df %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = condition), size=3, pch=21) +
  labs(fill = "", title = "Sample Type") +
  scale_fill_manual(values = condition_colors) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  )

pre_plots <- (p1 | p2)
pre_plots

p3 <- post_df %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = study), size=3, pch=21) +
  labs(fill = "", title = "Study") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  )
p4 <- post_df %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = condition), size=3, pch=21) +
  labs(fill = "", title = "Sample Type") +
  scale_fill_manual(values = condition_colors) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5)
  )

post_plots <- (p3 | p4)
post_plots

study_plots <- p1 | p3
condition_plots <- p2 | p4

# mainFigure <- (pre_plots / post_plots) + plot_layout(guides = "collect")
mainFigure <- (study_plots / condition_plots) + plot_layout(guides = "collect")


# Code for UMAP plot
library(matrixStats)
library(umap)

gene_vars <- rowVars(assay(post_vsd))
head(gene_vars)
names(gene_vars) <- rownames(post_vsd)

ngenes = 1000
hvgs <- sort(gene_vars) %>% tail(ngenes) %>% names()

pca_res <- prcomp(t(assay(post_vsd)[hvgs, ]))

set.seed(4132)
n_dims <- 50
umap_res <- umap(pca_res$x[, 1:n_dims])


umap_df <- data.frame(
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  condition = post_vsd$condition,
  study = post_vsd$study
)


umapPlot <- umap_df %>%
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(fill = condition), pch=21, size = 3) +
  labs(fill = "") +
  scale_fill_manual(values = condition_colors) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 14),
    legend.position = 'bottom'
  )
umapPlot


## Save for later usage:
coord_df <- pre_df %>%
  rename(pre.PC1 = PC1, pre.PC2 = PC2) %>%
  select(name, study, condition, pre.PC1, pre.PC2) %>%
  inner_join(
    post_df %>% 
      rename(post.PC1 = PC1, post.PC2 = PC2) %>% 
      select(name, study, condition, post.PC1, post.PC2),
    by = c("name", "condition"),
  ) %>%
  inner_join(
    umap_df %>% tibble::rownames_to_column("name") %>% select(-study),
    by = c("name", "condition")
  )
write_csv(coord_df, file.path("data", "manuscript", "Figure1_Coordinate_Data.csv"))

## Save figures
ggsave(
  filename = file.path(figureDir, "Figure1_PCA_Plots.svg"),
  plot = mainFigure,
  width = 10,
  height = 6
)
ggsave(
  filename = file.path(figureDir, "Figure1_UMAP_Plot.svg"),
  plot = umapPlot,
  width = 5,
  height = 6
)
