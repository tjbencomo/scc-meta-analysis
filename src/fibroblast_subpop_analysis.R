## Description: Analyze subpopulations of fibroblasts in Andrew Ji's scRNA-seq data
## using Schutz 2023 signatures

library(Seurat)
library(dplyr)
library(UCell)
library(readxl)
library(readr)
library(patchwork)
library(ggplot2)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
schutzDir <- file.path(dataDir, "schutz-2023")
cells <- readRDS(file.path(dataDir, "ji-2019", "Ji_2019_Cells.rds"))

cells@meta.data %>%
  count(level1_celltype) %>%
  View()
cells@meta.data %>%
  count(level2_celltype) %>%
  View()
cells@meta.data %>%
  count(level3_celltype) %>%
  View()


fibs <- subset(cells, level2_celltype == "Fibroblast")

fibs <- FindVariableFeatures(fibs)
fibs <- ScaleData(fibs)
fibs <- RunPCA(fibs)
fibs <- FindNeighbors(fibs)
fibs <- RunUMAP(fibs, dims = 1:50)

DimPlot(fibs, group.by = "tum.norm")
fibs <- FindClusters(fibs, resolution = .3)
umapClusterPlot <- DimPlot(fibs)
umapClusterPlot

proInflamSig <- read_excel(file.path(schutzDir, "normal_skin_fibroblast_populations.xlsx"), 
                           sheet = "Cluster 1 (pro-inflammatory)", skip = 2) %>% pull(gene)
secretReticSig <- read_excel(file.path(schutzDir, "normal_skin_fibroblast_populations.xlsx"), 
                             sheet = "Cluster 3 (secretory-reticular)") %>% pull(gene)
mesenchSig <- read_excel(file.path(schutzDir, "normal_skin_fibroblast_populations.xlsx"), 
                         sheet = "Cluster 6 (mesenchymal)") %>% pull(gene)
secretPapSig <- read_excel(file.path(schutzDir, "normal_skin_fibroblast_populations.xlsx"), 
                           sheet = "Cluster 13 (secretory-papillary") %>% pull(gene)
icafSig <- read_excel(file.path(schutzDir, "sccis_scc_fibroblast_populations.xlsx"), 
                      sheet = "iCAF", skip = 2) %>% pull(gene)
mycafSig <- read_excel(file.path(schutzDir, "sccis_scc_fibroblast_populations.xlsx"), 
                      sheet = "myCAF") %>% pull(gene)

gs <- list(
  'Pro_Inflammatory' = proInflamSig,
  'Secretory_Reticular'=  secretReticSig,
  'Mesenchymal' = mesenchSig,
  'Secretory_Papillary' = secretPapSig,
  'iCAF' = icafSig,
  'myCAF' = mycafSig
)

fibs <- AddModuleScore_UCell(fibs, features = gs)

VlnPlot(fibs, paste0(names(gs), "_UCell"))

scorePlots <- VlnPlot(fibs, c("iCAF_UCell", "myCAF_UCell"), pt.size = 0)
scorePlots

sigs <- paste0(names(gs), "_UCell")
VlnPlot(fibs, sigs, pt.size = 0)

degs <- FindAllMarkers(fibs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
topMarkers <- degs %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

degs %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  filter(cluster == 4) %>%
  pull(gene) %>%
  write_clip()

# 2: inflammatory iCAF + secretory reticular
# 0/1/3/4: mesenchymal + myCAF + secretory papillary
# Note Enrichr on top 100 marker genes for cluster 4 = smooth muscle contraction=myCAF

fib_annot_df <- fibs@meta.data %>%
  mutate(fibroblast_cluster_id = Idents(fibs)) %>%
  mutate(fibroblast_type = case_when(
    fibroblast_cluster_id == 2 ~ "Inflam Fibroblast ",
    fibroblast_cluster_id %in% c(0,1,3,4) ~ "Mesenchymal Fibroblast"
  )) %>%
  tibble::rownames_to_column("barcode")
write_csv(fib_annot_df, "data/ji-2019/fibroblast_annotations.csv")

stopifnot(all(fib_annot_df$barcode == colnames(fibs)))
fibs$label <- stringr::str_replace(fib_annot_df$fibroblast_type, " Fibroblast", "")
umapLabelPlot <- DimPlot(fibs, group.by = "label") + labs(title = "")
umapLabelPlot

umapPlots <- umapClusterPlot | umapLabelPlot
umapPlots


## Save plots
ggsave(
  filename = file.path(figureDir, "Figure4_Fibroblast_UMAPs.svg"),
  plot = umapPlots,
  width = 8,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_Fibroblast_Sig_Scores.svg"),
  plot = scorePlots,
  width = 6, height = 4
)
