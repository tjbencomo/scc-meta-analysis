## Description: Draw heatmap showing expression of 
## genes nominated by Shain meta analysis as driver genes

library(DESeq2)
library(readr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")

## Shain driver genes
# Genes pulled from paper
shainGenes <- c("TP53", "NOTCH1", "CDKN2A", "FAT1", "ARID2", "CASP8", "NOTCH2", 
                "HRAS", "TPO", "CHUK", "PBRM1", "USP28", "COL11A1", "CADPS", 
                "NOVA1", "SORCS3", "DMD", "DNAH7", "COL4A5", "CEP85L", 
                "DGKI", "DOCK9", "CREBBP", "AJUBA", "PTEN", "PIK3CA", "EZH2",
                "KRAS", "CCND1", "MTOR")
# From supplementary data table 2
# shainGenes <- c(
#   CHUK, PBRM1, USP28, COL11A1, CADPS, NOVA1, TPO, DMD, DNAH7, COL4A5,
#   CEP85L, HRAS, TP53, FAT1, ARID2, CDKN2A, NOTCH1, NOTCH2, CASP8,
#   SORCS3, WHSC1, DGKI, DOCK9
# )
resdf <- read_csv(file.path(deseqDir, "SCC_vs_NS.csv"))

vsd <- readRDS(file.path(deseqDir, "limma_batch_normalized.rds"))
ensgene2symbol <- rowData(vsd)$symbol
names(ensgene2symbol) <- rownames(rowData(vsd))

idx <- which(ensgene2symbol %in% shainGenes)
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
  summarise(across(all_of(shainGenes), list(mean))) %>%
  filter(condition != "AK_IEC", condition != "AK_IEC_SCC") %>%
  as.data.frame() %>%
  tibble::column_to_rownames("condition") %>%
  t()
rownames(avgMat) <- str_replace(rownames(avgMat), "_1", "")
 
condition_order <- c("NS", "AK", "IEC", "SCC", "KA")

avgMat2 <- avgMat %>%
  t() %>%
  scale() %>%
  t()


p <- Heatmap(
  avgMat2[, condition_order],
  cluster_columns = FALSE,
  show_row_dend = FALSE
)
p

svg(file.path(figureDir, "FigureS3_Driver_Gene_Expression.svg"), width = 5, height = 5)
p
dev.off()
