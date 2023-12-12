## Description: Perform differential expression analysis for SCC vs NS for each study
## individually. Save results and calculate gene-level statistics across studies
## Used to understand how DEGs vary across studies
## 

library(DESeq2)
library(IHW)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(viridis)


dataDir <- file.path("data")
deseqDir <- file.path("data", "processed", "deseq")

dds <- readRDS(file.path(deseqDir, "deseq_obj.rds"))
pooled_df <- read_csv(file.path(deseqDir, "SCC_vs_NS.csv"))

# Figure out which studies have NS > 1 and SCC > 1 for individual DEG analysis
x <- table(dds$study, dds$condition)
nsIdx <- x[, 1] > 1
sccIdx <- x[, 7] > 1
keepIdx <- nsIdx & sccIdx
studies <- rownames(x)[keepIdx]

dds <- dds[, dds$study %in% studies]

# studies <- studies[1:2]
df_list <- list()
for (i in 1:length(studies)) {
  study <- studies[i]
  print(paste("Analyzing", study))
  x <- dds[, dds$study == study & dds$condition %in% c("NS", "SCC")]
  x$condition <- factor(x$condition, levels = c("NS", "SCC"))
  x$study <- factor(x$study)
  design(x) <- ~ condition
  x <- DESeq(x)
  res <- results(x, contrast = c("condition", "SCC", "NS"), filterFun = ihw)
  resdf <- data.frame(res)
  resdf$gene <- rowData(x)$symbol
  resdf$ensgene <- rownames(resdf)
  df_list[[study]] <- resdf
}

resdf <- bind_rows(df_list, .id = "study") %>% as_tibble()
rm(dds, df_list, x, res); gc()

# Question: Should we filter out studies with log2FC of NA to get info on more
# genes?
stats_df <- resdf %>%
  group_by(ensgene, gene) %>%
  summarize(
    mean_expression = mean(baseMean),
    log2FC_max = max(log2FoldChange),
    log2FC_min = min(log2FoldChange),
    log2FC_mean = mean(log2FoldChange),
    log2FC_sd = sd(log2FoldChange),
    pval_max = max(pvalue),
    pval_min = min(pvalue),
    pval_mean = mean(pvalue),
    pval_sd = sd(pvalue),
    padj_max = max(padj),
    padj_min = min(padj),
    padj_mean = mean(padj),
    padj_sd = sd(padj)
  ) %>%
  mutate(both_directions = sign(log2FC_max) != sign(log2FC_min)) %>%
  ungroup()



# Save individual study DEG results
write_csv(resdf, file.path(dataDir, "study_level_deg_estimates.csv.gz"))
# Save gene-level statistics
write_csv(stats_df, file.path(dataDir, "study_level_gene_statistics.csv.gz"))

