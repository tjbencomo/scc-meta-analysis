## Description: Plot top DEGs upregulated in SCC compared to NS using scRNA data from
## Andrew's cSCC samples to show cell type specificity

library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(readxl)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path("data", "processed", "deseq")

resdf <- read_csv(file.path(deseqDir, "SCC_vs_NS.csv"))
priorGenesDf <- read_excel(file.path(dataDir, "Previous_cSCC_DEGs.xlsx", sheet = 2))
cells <- readRDS(file.path(dataDir, "ji-2019", "Ji_2019_Cells.rds"))


cells$collapsed_level2 <- with(
  cells@meta.data,
  case_when(
    level2_celltype %in% c("ASDC", "CD1C", "CLEC9A", "PDC") ~ "Dendritic Cell",
    level2_celltype %in% c("LC", "Mac") ~ "Macrophage",
    level2_celltype == "MDSC" ~ "MDSC",
    level2_celltype %in% c("ASDC", "CD1C", "CLEC9A", "LC") ~ "Dendritic Cell",
    TRUE ~ level2_celltype
  )
)

cells <- subset(cells, level2_celltype != "Keratinocyte")

tumor <- subset(cells, tum.norm == "Tumor")
normal <- subset(cells, tum.norm == "Normal")

tumor_level2_col_order <- c("Tumor_KC_Basal", "Tumor_KC_Cyc", "Tumor_KC_Diff", "TSK", 
                      "Pilosebaceous", "Eccrine",
                      "Fibroblast", "Endothelial Cell",
                      "Tcell", "B Cell", "NK",
                      "Dendritic Cell", "Macrophage", "MDSC",
                      "Melanocyte")
normal_level2_col_order <- c("Normal_KC_Basal", "Normal_KC_Cyc", "Normal_KC_Diff",
                            "Pilosebaceous", "Eccrine",
                            "Fibroblast", "Endothelial Cell",
                            "Tcell", "B Cell", "NK",
                            "Dendritic Cell", "Macrophage", "MDSC",
                            "Melanocyte")
tumor$collapsed_level2 <- factor(tumor$collapsed_level2, levels = tumor_level2_col_order)
normal$collapsed_level2 <- factor(normal$collapsed_level2, levels = normal_level2_col_order)


priorGenes <- priorGenesDf %>% 
  filter(Direction == "Up") %>%
  count(Gene) %>%
  filter(n > 1) %>%
  pull(Gene)
topBulkGenes <- resdf %>%
  slice_max(stat, n = 30)
upGeneTally <- read_csv(file.path(dataDir, "Up_Gene_Tally.csv"))

upGenes <- unique(c(priorGenes, topBulkGenes$gene_symbol, upGeneTally$gene_symbol[upGeneTally$n_studies == 9]))
upGeneOrder <- resdf %>%
  filter(gene_symbol %in% upGenes) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(gene_symbol)

upPlot <- DotPlot_scCustom(
  seurat_object = tumor, 
  features = rev(upGeneOrder),
  group.by = "collapsed_level2",
  x_lab_rotate = TRUE,
  flip_axes = TRUE,
  colors_use = colorRampPalette(c("blue", "red"))(10)
)
upPlot
VlnPlot(tumor, "KRT16", group.by = "collapsed_level2")

dnGeneTally <- read_csv(file.path(dataDir, "Dn_Gene_Tally.csv"))
bottomBulkGenes <- resdf %>%
  slice_min(stat, n = 30)
downGenes <- unique(c(bottomBulkGenes$gene_symbol, "CCL27", dnGeneTally$gene_symbol[dnGeneTally$n_studies == 9]))
downGeneOrder <- resdf %>%
  filter(gene_symbol %in% downGenes) %>%
  arrange(log2FoldChange) %>%
  pull(gene_symbol)

downPlot <- DotPlot_scCustom(
  seurat_object = normal, 
  features = rev(downGeneOrder),
  group.by = "collapsed_level2",
  x_lab_rotate = TRUE,
  flip_axes = TRUE,
  colors_use = colorRampPalette(c("blue", "red"))(10)
)
downPlot

## Save plots
ggsave(
  filename = file.path(figureDir, "Figure3_UpGenes_scRNA.svg"),
  plot = upPlot,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(figureDir, "Figure3_DownGenes_scRNA.svg"),
  plot = downPlot,
  width = 12,
  height = 8
)
