## Examine UMAP and PCA plots to see what samples are outliers using ggplotly
## See Problem_Samples google sheets for details on bad samples

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(matrixStats)
library(umap)

dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")
vsd <- readRDS(file.path(deseqDir, "limma_batch_normalized.rds"))

gene_vars <- rowVars(assay(vsd))
head(gene_vars)
names(gene_vars) <- rownames(vsd)

ngenes = 500
hvgs <- sort(gene_vars) %>% tail(ngenes) %>% names()

pca_res <- prcomp(t(assay(vsd)[hvgs, ]))
umap_res <- umap(pca_res$x)

dim_df <- data.frame(
  sample_id = vsd$sample_id,
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = vsd$condition,
  study = vsd$study
)

label_ids <- c(290, 261, 336, 252, 54, 236, 53, 57, 3, 19, 54, 16)
label_samples <- paste0("Sample_", unique(label_ids))

pcaPlot <- dim_df %>%
  mutate(text_label = ifelse(sample_id %in% label_samples, sample_id, "")) %>%
  ggplot(aes(PC1, PC2, label = text_label)) +
  geom_point(aes(fill = condition), pch=21, size = 3) +
  theme_bw() +
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0) +
  scale_fill_brewer(palette = "Accent2")
pcaPlot
# ggplotly(pcaPlot)


umapPlot <- dim_df %>%
  mutate(text_label = ifelse(sample_id %in% label_samples, sample_id, "")) %>%
  ggplot(aes(UMAP1, UMAP2, label = text_label)) +
  geom_point(aes(fill = condition), pch=21, size = 3) +
  theme_bw() +
  geom_text_repel(max.overlaps = Inf, min.segment.length = 0) +
  scale_fill_brewer(palette = "Accent2")
umapPlot

# ggplotly(umapPlot)
