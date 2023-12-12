## Description: Draw figures for GTEx sun exposure analysis

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(enrichR)
library(fgsea)
library(patchwork)
library(stringr)

theme_set(theme_classic())


data_dir <- file.path("data")
figureDir <- file.path("figures", "manuscript")
map_df <- read_csv(file.path(data_dir, "Sun_Exposed_vs_Protected_DE_MAP.csv"))
vsd <- readRDS(file.path(data_dir, "processed", "deseq","limma_batch_normalized.rds"))
score_df <- read_csv(file.path(data_dir, "gsea-results", "Hallmark_Reactome_DvP_GSVA_Scores.csv"))
scc_df <- read_csv(file.path(data_dir, "processed", "deseq", "SCC_vs_NS.csv"))

condition_colors <- readRDS(file.path(data_dir, "condition_color_palette.rds"))
blured_colors <- readRDS(file.path(data_dir, "up_down_direction_colors.rds"))


## Decision to filter based on expression:
## Without threshold no Enrichr hits
evp_df <- map_df %>%
  filter(gene %in% rowData(vsd)$symbol) %>%
  group_by(up_in_exposed) %>%
  filter(
    padj < 1e-6, 
    baseMean > 100
  ) %>%
  slice_max(abs(log2FoldChange), n = 150) %>%
  ungroup()
evpSig <- evp_df$log2FoldChange
names(evpSig) <- evp_df$gene

write_csv(evp_df, file.path(data_dir, "GTEx_EvP_Signature.csv"))

## Volcano Plot highlighting signature genes
showGenes <- c("MMP3", "OTOP3", "ZIC1", "KRTAP13-1", "DHRS2", "MAB21L1", "PRG4", "KRT38")

gtexMaPlot <- map_df %>%
  mutate(
    color_label = case_when(
      gene %in% evp_df$gene & log2FoldChange > 0 ~ "up",
      gene %in% evp_df$gene & log2FoldChange < 0 ~ "down",
      TRUE ~ "ignore"
    )
  ) %>%
  mutate(color_label = factor(color_label, levels=c("ignore", "down", "up"))) %>%
  arrange(color_label) %>%
  mutate(text_label = case_when(
    gene %in% showGenes ~ gene,
    TRUE ~ ""
  )) %>%
  ggplot(aes(log10(baseMean), log2FoldChange, label = text_label)) +
  geom_point(aes(fill = color_label), pch = 21, alpha = .7) +
  guides(fill = "none") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue", "ignore" = "gray")) +
  geom_hline(yintercept = c(-.5, .5), linetype = "dashed", color = "red") +
  geom_text_repel(min.segment.length = 0, max.overlaps = Inf, box.padding = 1) +
  scale_y_continuous(breaks = seq(-8, 5, by = 2)) +
  labs(x = "Log10 Mean Expression", y = "Log2 Fold Change") +
  theme(text = element_text(size = 14))
gtexMaPlot


## Enrichr on EvP signature
# setEnrichrSite("Enrichr")
# db <- c("MSigDB_Hallmark_2020")
# up_enrichr_df <- enrichr(evp_df$gene[evp_df$log2FoldChange > 0], db)
# dn_enrichr_df <- enrichr(evp_df$gene[evp_df$log2FoldChange < 0], db)
# 
# enrichr_df <- bind_rows(
#   up_enrichr_df$MSigDB_Hallmark_2020 %>%
#     mutate(direction = "up"),
#   dn_enrichr_df$MSigDB_Hallmark_2020 %>%
#     mutate(direction = "down")
# ) %>% as_tibble()
# 
# enrichPlotUp <- enrichr_df %>%
#   filter(Adjusted.P.value < .1, direction == "up") %>%
#   ggplot(aes(Combined.Score, reorder(Term, Combined.Score))) +
#   geom_col(fill = "red", alpha = .7) +
#   labs(x = "Combined Score", y = "", title = "Up in Sun Exposure") +
#   theme(plot.title = element_text(hjust = .5))
# enrichPlotUp
# enrichPlotDn <- enrichr_df %>%
#   filter(Adjusted.P.value < .1, direction == "down") %>%
#   ggplot(aes(Combined.Score, reorder(Term, Combined.Score))) +
#   geom_col(fill = "blue", alpha = .7) +
#   labs(x = "Combined Score", y = "") +
#   scale_y_discrete(position = "right") +
#   scale_x_reverse() +
#   labs(x = "Combined Score", y = "", title = "Down in Sun Exposure") +
#   theme(plot.title = element_text(hjust = .5))
# enrichPlotDn
# 
# enrichPlot <- enrichPlotDn | enrichPlotUp

## Score samples with EvP signature
rownames(vsd) <- rowData(vsd)$symbol
geneMat <- assay(vsd)[evp_df$gene, ]
evpScores <-  colSums(geneMat * evp_df$log2FoldChange)
vsd$EvP_Score <- scale(evpScores)[, 1]


combined_df <- colData(vsd) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(sampleID = str_c(sample_id, condition, sep="-")) %>%
  inner_join(
    score_df %>% 
      select(sampleID, DvP_Score, HALLMARK_UV_RESPONSE_UP, HALLMARK_UV_RESPONSE_DN)
  )
write_csv(combined_df, file.path(data_dir, "gsea-results", "EvP_Scores.csv"))


evpPlot <- combined_df %>%
  filter(condition != "AK_IEC", condition != "AK_IEC_SCC") %>%
  mutate(condition = factor(condition, levels=c("NS", "AK", "IEC", "SCC", "KA"))) %>%
  ggplot(aes(condition, EvP_Score)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "EvP Score", fill = "") +
  theme(text = element_text(size = 14)) +
  guides(fill = "none")
evpPlot
lm(EvP_Score ~ condition, data = colData(vsd)) %>% summary()
kruskal.test(EvP_Score ~ condition, data = colData(vsd))


## Compare EvP in RDEB vs sporadic cSCC
rdeb_df <- combined_df %>%
  mutate(RDEB_Status = case_when(
    condition == "SCC" & study == "Cho_2018" ~ "RDEB",
    condition == "SCC" & study != "Cho_2018" ~ "Sporadic",
    TRUE ~ "Other"
  )) %>%
  filter(condition == "SCC")
rdebPlot <- rdeb_df %>%
  ggplot(aes(RDEB_Status, EvP_Score)) +
  geom_boxplot(aes(fill = RDEB_Status)) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = "", y = "EvP Score") +
  guides(fill = "none") +
  theme(text = element_text(size = 14))
rdebPlot
wilcox.test(EvP_Score ~ RDEB_Status, data = rdeb_df)


## Show EvP compared to SCC vs NS genes
## What DE genes show same direction and what show different changes
conc_disc_colors <- blured_colors
names(conc_disc_colors) <- c("Discordant", "Concordant")
conc_disc_colors <- c(conc_disc_colors, c("ns" = "gray"))
comparisonPlot <- evp_df %>%
  inner_join(
    scc_df, by = c("gene" = "gene_symbol"),
    suffix = c(".gtex", ".scc")
  ) %>%
  mutate(
    color_label = case_when(
      sign(log2FoldChange.gtex) == sign(log2FoldChange.scc) & padj.scc < .05 ~ "Concordant",
      sign(log2FoldChange.gtex) != sign(log2FoldChange.scc) & padj.scc < .05 ~ "Discordant",
      TRUE ~ "ns"
    )
  ) %>%
  group_by(log2FoldChange.gtex > 0, log2FoldChange.scc > 0) %>%
  arrange(desc(abs(log2FoldChange.gtex + log2FoldChange.scc))) %>%
  mutate(rank_id = row_number()) %>%
  ungroup() %>%
  mutate(
    text_label = case_when(
      rank_id < 6 & color_label != "ns" ~ gene,
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(log2FoldChange.gtex, log2FoldChange.scc, label = text_label)) +
  geom_point(aes(fill = color_label), pch=21, size=2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = conc_disc_colors) +
  # scale_fill_manual(values = c("Concordant" = "green", "Discordant" = "red", "ns" = "gray")) +
  geom_text_repel(min.segment.length = 0, max.overlaps = Inf) +
  labs(x = "Sun Exposed vs Protected Log2FC", y = "SCC vs NS Log2FC", fill = "") +
  theme(legend.position = 'bottom',
        text = element_text(size = 14))
comparisonPlot

p1 <- (gtexMaPlot | evpPlot) + plot_layout(widths = c(2, 1))
p2 <- (rdebPlot | comparisonPlot) + plot_layout(widths = c(1, 3))

finalPlot <- p1 / p2



## Save figures
ggsave(
  filename = file.path(figureDir, "Figure5_Sun_Exposed_vs_Protected_MA.png"),
  plot = gtexMaPlot,
  width = 5,
  height = 5
)


ggsave(
  filename = file.path(figureDir, "Figure5_EvP_Condition_Scores.svg"),
  plot = evpPlot,
  width = 5,
  height = 4
)


ggsave(
  filename = file.path(figureDir, "Figure5_RDEB_EvP_Scores.svg"),
  plot = rdebPlot,
  width = 3,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure5_DE_Comparisons.svg"),
  plot = comparisonPlot,
  width = 6,
  height = 6
)

ggsave(
  filename = file.path(figureDir, "Figure5_All_Panels.svg"),
  plot = finalPlot,
  width = 7.5,
  height = 10
)



