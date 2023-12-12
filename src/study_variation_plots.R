## Description: Plot figures for study-level DEG variation analysis

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(patchwork)
library(UpSetR)
library(ggrepel)
library(tidyr)
library(enrichR)
library(stringr)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path('data')
deseqDir <- file.path("data", "processed", "deseq")
stats_df <- read_csv(file.path(dataDir, 'study_level_gene_statistics.csv.gz'))
study_df <- read_csv(file.path(dataDir, 'study_level_deg_estimates.csv.gz'))
pooled_df <- read_csv(file.path(dataDir, 'processed', 'deseq', 'SCC_vs_NS.csv'))
met_df <- read_csv(file.path(dataDir, "metadata_final_cohort.csv"))
vsd <- readRDS(file.path(deseqDir, "limma_batch_normalized.rds"))
up_dn_colors <- readRDS(file.path(dataDir, "up_down_direction_colors.rds"))

geneInfo <- data.frame(gene_id = rownames(vsd), gene_symbol = rowData(vsd)$symbol)

# Compare DEGs between studies using FDR < 5%
deg_list <- list()
up_deg_list <- list()
dn_deg_list <- list()
studies <- unique(study_df$study)
for (i in 1:length(studies)) {
  study_name <- studies[i]
  print(paste("Getting DEGs for", study_name))
  deg_list[[study_name]] <- study_df %>%
    filter(study == study_name, padj < .05) %>%
    pull(ensgene)
  up_deg_list[[study_name]] <- study_df %>%
    filter(study == study_name, padj < .05, log2FoldChange > 0) %>%
    pull(ensgene)
  dn_deg_list[[study_name]] <- study_df %>%
    filter(study == study_name, padj < .05, log2FoldChange < 0) %>%
    pull(ensgene)
}
names(up_deg_list) <- str_replace(names(up_deg_list), '_', ' ')
names(dn_deg_list) <- str_replace(names(dn_deg_list), '_', ' ')

upOverlapPlot <- upset(fromList(up_deg_list), order.by = 'degree', nsets = length(up_deg_list), 
                       nintersects = 35, text.scale = 2)
print(upOverlapPlot)

dnOverlapPlot <- upset(fromList(dn_deg_list), order.by = 'degree', nsets = length(dn_deg_list),
                       nintersects = 35, text.scale = 2)
print(dnOverlapPlot)

up_gene_tally_df <- data.frame(
  n_studies = table(unlist(up_deg_list))
) %>% rename(n_studies = n_studies.Freq, gene_id = n_studies.Var1) %>%
  arrange(desc(n_studies)) %>%
  inner_join(geneInfo)
dn_gene_tally_df <- data.frame(
  n_studies = table(unlist(dn_deg_list))
) %>% rename(n_studies = n_studies.Freq, gene_id = n_studies.Var1) %>%
  arrange(desc(n_studies)) %>%
  inner_join(geneInfo)

# This table tallies number of studies each gene is up/down in
# For specific studies see study_level_deg_estimates.csv.gz
tally_df <- up_gene_tally_df %>%
  rename(n_upregulated_studies = n_studies) %>%
  full_join(
    dn_gene_tally_df %>% rename(n_downregulated_studies = n_studies)
  ) %>%
  replace_na(list(n_upregulated_studies = 0, n_downregulated_studies = 0))
write_csv(tally_df, "data/gene_de_study_counts.csv.gz")

## Save these results
write_csv(up_gene_tally_df, file.path(dataDir, "Up_Gene_Tally.csv"))
write_csv(dn_gene_tally_df, file.path(dataDir, "Dn_Gene_Tally.csv"))

studyCounts <- sapply(deg_list, length)
counts_df <- data.frame(
  study = names(studyCounts),
  n_deg = studyCounts
) %>% as_tibble()

degCountPlot <- counts_df %>%
  mutate(study = case_when(
    study == "Mahapatra_2020" ~ "Das Mahapatra_2020",
    TRUE ~ study
  )) %>%
  mutate(study = str_replace_all(study, "_", "\n")) %>%
  ggplot(aes(n_deg, reorder(study, n_deg))) +
  geom_col() +
  theme_classic() +
  # geom_text(aes(label = n_deg), hjust = -0.1) +
  labs(x = "Number of DEGs", y = "") +
  theme(plot.title = element_text(hjust = .5),
        text = element_text(size = 14)) +
  scale_x_continuous(breaks = seq(0, 18000, by = 2000), limits = c(0, 18000))
degCountPlot

study_size_deg_df <- counts_df %>%
  inner_join(
    met_df %>%
      count(study_name, name = "n_samples") %>%
      select(study_name, n_samples),
    by = c("study" = "study_name")
  )
studySizePlot <- study_size_deg_df %>%
  mutate(study = str_replace(study, "_", " ")) %>%
  ggplot(aes(n_samples, n_deg, label = study)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_smooth(method = 'lm') +
  labs(x = "Sample Size", y = "Number of DEGs") +
  scale_y_log10() +
  scale_x_log10() +
  geom_text_repel() +
  labs(x = "Log10 [Sample Size]", y = "Log10 [Number of DEGs]") +
  theme(
    text = element_text(size = 14)
  )
studySizePlot
cor.test(log10(study_size_deg_df$n_samples), log10(study_size_deg_df$n_deg), method = 'pearson')

# Repeat comparison using only canonical genes
# Compare DEGs between studies using FDR < 5%
up_genes <- c("CDKN2A", "FN1", "KRT16", "KRT17", "MMP1", "MMP10", "PI3",
              "PTHLH", "S100A12")
down_genes <- c("CCL27")
canon_deg_list <- list()
studies <- unique(study_df$study)
for (i in 1:length(studies)) {
  study_name <- studies[i]
  print(paste("Getting DEGs for", study_name))
  canon_deg_list[[study_name]] <- study_df %>%
    filter(study == study_name, padj < .05, gene %in% c(up_genes, down_genes)) %>%
    pull(ensgene)
}

upset(fromList(canon_deg_list), order.by = 'degree', nsets = length(deg_list), text.scale = 2)
upset(fromList(canon_deg_list), order.by = 'freq', nsets = length(deg_list), text.scale = 2)

gene_hits <- data.frame(
  gene_id = unlist(canon_deg_list)
  ) %>%
  count(gene_id, name = "n_studies") %>%
  arrange(desc(n_studies)) %>%
  inner_join(
    select(pooled_df, gene_id, gene_symbol)
  )

gene_hits %>%
  select(gene_symbol, n_studies)

## Show p-value distributions for canonical genes
canon_genes <- c(up_genes, down_genes)
canonPvalPlot <- study_df %>%
  filter(gene %in% canon_genes) %>%
  ggplot(aes(gene, padj)) +
  geom_boxplot(aes(fill = gene), outlier.shape = NA) +
  geom_point() +
  geom_hline(yintercept = .05, linetype = "dashed", color = "red") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
    text = element_text(size = 14)
  ) +
  guides(fill = "none") +
  labs(x = "", y = "Adjusted P-Value") +
  scale_y_continuous(breaks = seq(0, .8, by = .1), limits = c(0, .8))
canonPvalPlot

## Draw plot comparing with microarray meta-analysis from Van Haren
microUpGenes <- c("VEGFC", "FN1", "WNT5A", "YY1", "MMP1", "CDH3")
microDnGenes <- c("KRT15")

microMetaPlot <- up_gene_tally_df %>%
  filter(gene_symbol %in% c(microUpGenes, microDnGenes)) %>%
  select(gene_symbol, n_studies) %>%
  rename(up_studies = n_studies) %>%
  full_join(
    dn_gene_tally_df %>%
      filter(gene_symbol %in% c(microUpGenes, microDnGenes)) %>%
      select(gene_symbol, n_studies) %>%
      rename(dn_studies = n_studies)
  ) %>%
  replace_na(list(up_studies = 0, dn_studies = 0)) %>%
  pivot_longer(!gene_symbol, names_to = "direction", values_to = "n_studies") %>%
  ggplot(aes(reorder(gene_symbol, -n_studies), n_studies, fill = direction)) +
  geom_bar(stat="identity", width=.75, position = "dodge", ) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  labs(x = "Gene", y = "Number of Studies", fill = "") +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) +
  scale_fill_manual(
    labels = c("dn_studies" = "Down in cSCC", "up_studies" = "Up in cSCC"),
    values = c("dn_studies" = up_dn_colors[[2]], "up_studies" = up_dn_colors[[1]])
  )
microMetaPlot

# pthlh_df <- data.frame(
#   PTHLH = assay(vsd)["ENSG00000087494.16", ],
#   study = vsd$study,
#   condition = vsd$condition
# )
mmp1_df <- data.frame(
  MMP1 = assay(vsd)["ENSG00000196611.6", ],
  study = vsd$study,
  condition = vsd$condition
)

condition_colors <- readRDS("data/condition_color_palette.rds")
mmp1Plot <- mmp1_df %>%
  filter(condition %in% c("NS", "SCC"), study %in% study_df$study) %>%
  mutate(study = case_when(
    study == "Mahapatra_2020" ~ "Das_Mahapatra_2020",
    TRUE ~ study
  )) %>%
  mutate(study = str_replace_all(study, "_", " ")) %>%
  ggplot(aes(condition, MMP1)) +
  geom_boxplot(aes(fill = condition)) +
  facet_wrap(~study, nrow = 3) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  labs(x = "", y = "Log2 MMP1", fill = "") +
  guides(fill = "none") +
  scale_fill_manual(values = condition_colors)
mmp1Plot

## Do EnrichR analysis
setEnrichrSite("Enrichr")
upGenes <- up_gene_tally_df %>%
  filter(n_studies > 6) %>%
  pull(gene_symbol)
dnGenes <- dn_gene_tally_df %>%
  filter(n_studies > 6) %>%
  pull(gene_symbol)
db <- c("MSigDB_Hallmark_2020")
up_enrichr_df <- enrichr(upGenes, db)
dn_enrichr_df <- enrichr(dnGenes, db)

enrichr_df <- bind_rows(
  up_enrichr_df$MSigDB_Hallmark_2020 %>%
    mutate(direction = "up"),
  dn_enrichr_df$MSigDB_Hallmark_2020 %>%
    mutate(direction = "down")
) %>% as_tibble()

enrichrPlot <- enrichr_df %>%
  filter(Adjusted.P.value < .1, Term != "Coagulation") %>%
  mutate(score = case_when(
    direction == "up" ~ -log10(Adjusted.P.value),
    direction == "down" ~ log10(Adjusted.P.value)
  )) %>%
  ggplot(aes(score, reorder(Term, score))) +
  geom_col(aes(fill = direction)) +
  theme_classic() +
  scale_fill_manual(
    labels = c("down" = "Down in cSCC", "up" = "Up in cSCC"),
    values = up_dn_colors
  ) +
  labs(x = "Log10 [Adjusted P-Value]", y = "", fill = "") +
  scale_x_continuous(limits = c(-4, 10), breaks = seq(-4, 10, by = 2))
enrichrPlot


## Save figures
svg(file.path(figureDir, "Figure2_Shared_Up_Genes.svg"), width = 18, height = 8)
upOverlapPlot
dev.off()

svg(file.path(figureDir, "Figure2_Shared_Dn_Genes.svg"), width = 14, height = 8)
dnOverlapPlot
dev.off()

ggsave(
  filename = file.path(figureDir, "Figure2_Number_DEGs.svg"),
  plot = degCountPlot,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure2_Mean_Variance.svg"),
  plot = expressionSdPlot,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure2_Sample_Size_DEG_Scatter.svg"),
  plot = studySizePlot,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure2_Canonical_Gene_Pvals.svg"),
  plot = canonPvalPlot,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure2_Microarray_MetaAnalysis_Comparison.svg"),
  plot = microMetaPlot,
  width = 7, height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure2_MMP1_Study_Variation.svg"),
  plot = mmp1Plot,
  width = 6, height = 6
)
ggsave(
  filename = file.path(figureDir, "Figure2_Enrichr_Top_Hits.svg"),
  plot = enrichrPlot,
  width = 6, height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure2_Top_Plots.svg"),
  plot = degCountPlot | canonPvalPlot | mmp1Plot,
  width = 18,
  height = 6
)
