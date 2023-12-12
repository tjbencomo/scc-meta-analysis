## Description: GSEA and GSVA analysis on NS vs AK and AK vs cSCC
## GSEA results from gsea_analysis.R

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

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")
gseaDir <- file.path(dataDir, "gsea-results")


## GSEA Plots
ns_ak_hallmark <- read_csv(file.path(gseaDir, "NS_vs_AK_Hallmark_GSEA.csv"))
ns_ak_reactome <- read_csv(file.path(gseaDir, "NS_vs_AK_Reactome_GSEA.csv"))
ak_scc_hallmark <- read_csv(file.path(gseaDir, "AK_vs_SCC_Hallmark_GSEA.csv"))
ak_scc_reactome <- read_csv(file.path(gseaDir, "AK_vs_SCC_Reactome_GSEA.csv"))

all_df <- bind_rows(
  ns_ak_hallmark %>% mutate(db = "Hallmarks", comparison = "NS vs AK"),
  ns_ak_reactome %>% mutate(db = "Reactome", comparison = "NS vs AK"),
  ak_scc_hallmark %>% mutate(db = "Hallmarks", comparison = "AK vs SCC"),
  ak_scc_reactome %>% mutate(db = "Reactome", comparison = "AK vs SCC"),
)


ns_ak_label_pathways <- c(
  "REACTOME_INNATE_IMMUNE_SYSTEM",
  "REACTOME_CELL_CYCLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_UV_RESPONSE_DN",
  "REACTOME_KERATINIZATION"
)
## NS vs AK Plot
p1 <- all_df %>%
  filter(comparison == "NS vs AK") %>%
  mutate(plabel = ifelse(pathway %in% ns_ak_label_pathways, pathway, "")) %>%
  mutate(plabel = str_replace(plabel, "HALLMARK_|REACTOME_", "")) %>%
  # mutate(plabel = pathway) %>%
  mutate(color_grp = case_when(
    abs(NES) > 1 & -log10(padj) > -log10(.05) ~ db,
    TRUE ~ "NS"
  )) %>%
  ggplot(aes(NES, -log10(padj), label = plabel)) +
  geom_point(aes(color = color_grp)) +
  theme_classic() +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_text_repel(min.segment.length = 0, max.overlaps = Inf, size=3) +
  labs(x = "Normalized Enrichment Score", y = "Log10 Adjusted P-Value", 
       title = "NS vs AK", color = "Database") +
  theme(plot.title = element_text(hjust = .5, face = "bold"),
        text = element_text(size = 14)) +
  scale_color_manual(values=c("Hallmarks" = "#F8766D", "Reactome" = "#00BA38", "NS" = "grey"))

# ggplotly(p)


ak_scc_label_pathways <- c(
  "REACTOME_METABOLISM_OF_LIPIDS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_KERATINIZATION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_CELL_CYCLE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_KRAS_SIGNALING_UP"
  
)
p2 <- all_df %>%
  filter(comparison == "AK vs SCC") %>%
  mutate(plabel = ifelse(pathway %in% ak_scc_label_pathways, pathway, "")) %>%
  mutate(plabel = str_replace(plabel, "HALLMARK_|REACTOME_", "")) %>%
  # mutate(plabel = pathway) %>%
  mutate(color_grp = case_when(
    abs(NES) > 1 & -log10(padj) > -log10(.05) ~ db,
    TRUE ~ "NS"
  )) %>%
  ggplot(aes(NES, -log10(padj), label = plabel)) +
  geom_point(aes(color = color_grp)) +
  theme_classic() +
  geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_text_repel(min.segment.length = 0, max.overlaps = Inf, size=3) +
  labs(x = "Normalized Enrichment Score", y = "Log10 Adjusted P-Value", 
       title = "AK vs SCC", color = "Database") +
  theme(plot.title = element_text(hjust = .5, face = "bold"),
        text = element_text(size = 14)) +
  scale_color_manual(values=c("Hallmarks" = "#F8766D", "Reactome" = "#00BA38", "NS" = "grey"))
# ggplotly(p2)

gseaVolcanoPlot <- (p1 | p2) + plot_layout(guides = "collect")



## Save plots
ggsave(
  filename = file.path(figureDir, "Figure3_NS_to_SCC_GSEA_Volcanos.svg"),
  plot = gseaVolcanoPlot,
  width = 9,
  height = 4
)
