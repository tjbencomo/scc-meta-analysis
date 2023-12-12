## Description: Calculate sample-level Hallmark and Reactome pathway scores with GSVA

library(DESeq2)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)
library(plotly)
library(viridis)
library(GSVA)
library(readxl)


figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")
gseaDir <- file.path(dataDir, "gsea-results")

print("Computing GSVA scores")
vsd <- readRDS(file.path(deseqDir, "limma_batch_normalized.rds"))
rownames(vsd) <- rowData(vsd)$symbol

hallmark <- gmtPathways(file.path(dataDir, "genesets", "h.all.v7.4.symbols.gmt"))
reactome <- gmtPathways(file.path(dataDir, "genesets", "c2.cp.reactome.v7.4.symbols.gmt"))

early_df_sig <- read_excel(file.path(dataDir, "Bailey_2023_Signatures.xlsx"), sheet = 1)$GeneID
late_df_sig <- read_excel(file.path(dataDir, "Bailey_2023_Signatures.xlsx"), sheet = 2)$GeneID
prog_sig <- read_excel(file.path(dataDir, "Bailey_2023_Signatures.xlsx"), sheet = 3)$Gene_ID
dvp_gs <- list('Early_Diff' = early_df_sig, 'Late_Diff' = late_df_sig, 'Progenitor' = prog_sig)

gs <- c(dvp_gs, hallmark, reactome)

# GSVA Analysis
esMat <- gsva(assay(vsd), gs, verbose = T)
esdf <- data.frame(t(esMat)) %>% 
    tibble::rownames_to_column("sampleID") %>% 
    as_tibble() %>%
    inner_join(data.frame(sampleID = colnames(vsd), condition = vsd$condition))


## Score with DvP signature via GeneFu
dvp_df <- read_excel(file.path(dataDir, "Bailey_2023_Signatures.xlsx"), sheet = "DvP_Signature") %>%
    filter(Gene_symbol %in% rownames(vsd))

dvpMat <- assay(vsd)[dvp_df$Gene_symbol, ]
dvpScores <-  colSums(dvpMat * dvp_df$Coefficient)
stopifnot(all(names(dvpScores) == esdf$sampleID))
esdf$DvP_Score <- scale(dvpScores)[, 1]

write_csv(esdf, file.path(gseaDir, "Hallmark_Reactome_DvP_GSVA_Scores.csv"))