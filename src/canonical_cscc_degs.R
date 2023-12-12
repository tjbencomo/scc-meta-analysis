## Description: Draw heatmap showing expression of 
## genes previously identified as DE in cSCC vs NS in microarray studies

library(DESeq2)
library(readr)
library(dplyr)
library(pheatmap)
library(ggplot2)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")

## Previously identified cSCC vs NS genes from microarray studies
up_genes <- c("CDKN2A", "FN1", "KRT16", "KRT17", "MMP1", "MMP10", "PI3",
              "PTHLH", "S100A12")
down_genes <- c("CCL27")
showGenes <- c(up_genes, down_genes)

vsd <- readRDS(file.path(deseqDir, "limma_batch_normalized.rds"))
ensgene2symbol <- rowData(vsd)$symbol
names(ensgene2symbol) <- rownames(rowData(vsd))

idx <- which(ensgene2symbol %in% showGenes)
showGenesInfo <- ensgene2symbol[idx]

rnadf <- data.frame(
  sampleID = colnames(vsd),
  condition = vsd$condition,
  t(assay(vsd)[names(showGenesInfo), ])
) %>% as_tibble()
stopifnot(colnames(rnadf)[3:ncol(rnadf)] == names(showGenesInfo))
colnames(rnadf)[3:ncol(rnadf)] <- showGenesInfo

## Average heatmap
avgMat <- rnadf %>%
  group_by(condition) %>%
  summarize(
    PTHLH = mean(PTHLH),
    FN1 = mean(FN1),
    PI3 = mean(PI3),
    KRT17 = mean(KRT17),
    CDKN2A = mean(CDKN2A),
    S100A12 = mean(S100A12),
    MMP10 = mean(MMP10),
    KRT16 = mean(KRT16),
    MMP1 = mean(MMP1),
    CCL27 = mean(CCL27)
  ) %>% 
  filter(condition != "AK_IEC", condition != "AK_IEC_SCC") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("condition") %>%
  t()
condition_order <- c("NS", "AK", "IEC", "SCC", "KA")
p <- pheatmap(
  avgMat[, condition_order],
  scale = 'row',
  cluster_cols = F,
  treeheight_row = 0, 
  treeheight_col = 0
)

svg(file.path(figureDir, "Figure3_Canonical_Gene_Heatmap.svg"), width = 4, height = 4)
p
dev.off()
