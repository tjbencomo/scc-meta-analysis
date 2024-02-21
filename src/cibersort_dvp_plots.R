## Description: Draw plots showing correlation between DvP score, KA vs SCC, and
## CIBERSORT deconvolution scores


library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
nonimmunedf <- read_csv(file.path(dataDir, "cibersort", "Limma_NonImmune_SCC_KA_Deconv.csv"))
immunedf <- read_csv(file.path(dataDir, "cibersort", "Limma_LM22_SCC_KA_Deconv.csv"))
metdf <- read_csv(file.path(dataDir, "metadata_final_cohort.csv"))
gsvadf <- read_csv(file.path(dataDir, "gsea-results", "Hallmark_Reactome_DvP_GSVA_Scores.csv"))
condition_colors <- readRDS(file.path(dataDir, "condition_color_palette.rds"))

df <- nonimmunedf %>%
  select(-`P-value`, -Correlation, -RMSE, -`Absolute score (sig.score)`) %>%
  inner_join(immunedf) %>%
  inner_join(metdf, by = c("Mixture" = "sample_id")) %>%
  inner_join(
    gsvadf %>% 
      select(sampleID, DvP_Score) %>%
      mutate(Mixture = str_split(sampleID, "-", simplify = T)[, 1])
  )
df$condition <- factor(df$condition, levels = c("SCC", "KA"))

theme_set(theme_classic())

## Compare cell types in KA vs SCC
basalCondPlot <- df %>%
  ggplot(aes(condition, Tumor_KC_Basal)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "Basal Tumor Cells") +
  guides(fill = "none")
basalCondPlot
wilcox.test(Tumor_KC_Basal ~ condition, data=df)

cycCondPlot <- df %>%
  ggplot(aes(condition, Tumor_KC_Cyc)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "Cycling Tumor Cells") +
  guides(fill = "none")
cycCondPlot
wilcox.test(Tumor_KC_Cyc ~ condition, data=df)

diffCondPlot <- df %>%
  ggplot(aes(condition, Tumor_KC_Diff)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "Differentiated Tumor Cells") +
  guides(fill = "none") +
  theme(
    text = element_text(size = 14)
  )
diffCondPlot
wilcox.test(Tumor_KC_Diff ~ condition, data=df)

tskCondPlot <- df %>%
  ggplot(aes(condition, TSK)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "TSKs") +
  guides(fill = "none") +
  theme(
    text = element_text(size = 14)
  )
tskCondPlot
wilcox.test(TSK ~ condition, data=df)

tumorSubpopCompPlots <- diffCondPlot | tskCondPlot
tumorSubpopCompPlots
# basalCondPlot | cycCondPlot | diffCondPlot | tskCondPlot

inflamFibCondPlot <- df %>%
  ggplot(aes(condition, Inflam_Fibroblast)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "Inflammatory Fibroblasts") +
  guides(fill = "none") +
  theme(
    text = element_text(size = 14)
  )
inflamFibCondPlot
wilcox.test(Inflam_Fibroblast ~ condition, data=df)

mesenchFibCondPlot <- df %>%
  ggplot(aes(condition, Mesenchymal_Fibroblast)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "", y = "Mesenchymal Fibroblasts") +
  guides(fill = "none") +
  theme(
    text = element_text(size = 14)
  )
mesenchFibCondPlot
wilcox.test(Mesenchymal_Fibroblast ~ condition, data=df)

fibSubpopCompPlots <- inflamFibCondPlot | mesenchFibCondPlot
fibSubpopCompPlots

## Compare tumor cell types in RDEB vs Sporadic cSCC
df <- df %>%
  mutate(RDEB_Status = case_when(
    condition == "SCC" & study_name == "Cho_2018" ~ "RDEB",
    condition == "SCC" & study_name != "Cho_2018" ~ "Sporadic",
    TRUE ~ "Other"
  ))

tskRdebPlot <- df %>%
  filter(condition == "SCC") %>%
  ggplot(aes(RDEB_Status, TSK)) +
  geom_boxplot(aes(fill = RDEB_Status)) +
  labs(x = "", y = "TSK") +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    text = element_text(size = 14)
  )
tskRdebPlot
wilcox.test(TSK ~ RDEB_Status, data=df %>% filter(condition == "SCC"))

diffRdebPlot <- df %>%
  filter(condition == "SCC") %>%
  ggplot(aes(RDEB_Status, Tumor_KC_Diff)) +
  geom_boxplot(aes(fill = RDEB_Status)) +
  labs(x = "", y = "Differentiated Tumor Cells") +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    text = element_text(size = 14)
  )
diffRdebPlot
wilcox.test(Tumor_KC_Diff ~ RDEB_Status, data=df %>% filter(condition == "SCC"))

mesenchRdebPlot <- df %>%
  filter(condition == "SCC") %>%
  ggplot(aes(RDEB_Status, Mesenchymal_Fibroblast)) +
  geom_boxplot(aes(fill = RDEB_Status)) +
  labs(x = "", y = "Mesenchymal Fibroblasts") +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    text = element_text(size = 14)
  )
mesenchRdebPlot
wilcox.test(Mesenchymal_Fibroblast ~ RDEB_Status, data=df %>% filter(condition == "SCC"))


rdebCompPlots <- diffRdebPlot | tskRdebPlot
rdebCompPlots

## Correlate DvP with cell types
dvpTskPlot <- df %>%
  ggplot(aes(DvP_Score, TSK)) +
  geom_point(size = 2) +
  # geom_smooth(method = 'lm') +
  labs(x = "DvP Score", y = "TSKs")
dvpTskPlot
cor.test(df$DvP_Score, df$TSK, method = 'spearman')

dvpDiffPlot <- df %>%
  ggplot(aes(DvP_Score, Tumor_KC_Diff)) +
  geom_point(size = 2) +
  # geom_smooth(method = 'lm') +
  labs(x = "DvP Score", y = "Differentiated Tumor Cells")
dvpDiffPlot
cor.test(df$DvP_Score, df$Tumor_KC_Diff, method = 'spearman')

dvpBasalPlot <- df %>%
  ggplot(aes(DvP_Score, Tumor_KC_Basal)) +
  geom_point(size = 2) +
  # geom_smooth(method = 'lm') +
  labs(x = "DvP Score", y = "Basal Tumor Cells")
dvpBasalPlot
cor.test(df$DvP_Score, df$Tumor_KC_Basal, method = 'spearman')

dvpCycPlot <- df %>%
  ggplot(aes(DvP_Score, Tumor_KC_Cyc)) +
  geom_point(size = 2) +
  # geom_smooth(method = 'lm') +
  labs(x = "DvP Score", y = "Cycling Tumor Cells")
dvpCycPlot
cor.test(df$DvP_Score, df$Tumor_KC_Cyc, method = 'spearman')


dvpInflamPlot <- df %>%
  ggplot(aes(DvP_Score, Inflam_Fibroblast)) +
  geom_point(size = 2) +
  labs(x = "DvP Score", y = "Inflammatory Fibroblasts")
dvpInflamPlot
cor.test(df$DvP_Score, df$Inflam_Fibroblast)

tskInflamPlot <- df %>%
  ggplot(aes(TSK, Inflam_Fibroblast)) +
  geom_point(size = 2) +
  labs(x = "TSK", y = "Inflammatory Fibroblasts")
tskInflamPlot
cor.test(df$TSK, df$Inflam_Fibroblast)

tskMesenchPlot <- df %>%
  ggplot(aes(TSK, Mesenchymal_Fibroblast)) +
  geom_point(size = 2) +
  labs(x = "TSK", y = "Mesenchymal Fibroblasts")
tskMesenchPlot
cor.test(df$TSK, df$Mesenchymal_Fibroblast)

infoCols <- c("Mixture", "P-value", "Correlation", "RMSE", "Absolute score (sig.score)")
cellTypeCols <- setdiff(c(colnames(nonimmunedf), colnames(immunedf)), infoCols)
deconMat <- as.matrix(df[df$condition == "SCC", cellTypeCols])

# ctCorPlot <- corrplot(cor(deconMat, method = 'spearman'), type = 'upper')

diffCd8Plot <- df %>%
  filter(condition == "SCC") %>%
  ggplot(aes(Tumor_KC_Diff, `T cells CD8`)) +
  geom_point(size = 2) +
  labs(x = "Differentiated Tumor Cells", y = "CD8 T Cells")
diffCd8Plot
with(
  df %>% filter(condition == "SCC"),
  cor.test(Tumor_KC_Diff, `T cells CD8`)
)

tskCd8Plot <- df %>%
  filter(condition == "SCC") %>%
  ggplot(aes(TSK, `T cells CD8`)) +
  geom_point(size = 2) +
  labs(x = "TSK", y = "CD8 T Cells")
tskCd8Plot
with(
  df %>% filter(condition == "SCC"),
  cor.test(TSK, `T cells CD8`)
)


## Compare between immunosuppressed (IS) and immunocompetent (IC)
## Veenstra_2023 and Chitsazzadeh_2016 excluded immunosuppressed patients = IC
mahapdf <- read_csv("data/Mahapatra_Immune_Info.csv")
baileydf <- read_csv("data/Bailey_Supplemental_File1.csv")
baileydf$original_sample_id <- str_c("PRNJA844527", baileydf$SEQ_ID, sep="_")
lukowskidf <- read_csv("data/E-MTAB-6430.csv")
lukowskidf <- lukowskidf %>% 
  rename(original_sample_id = sample_id, sample_accession = ena_accession)
  

immunodf <- df %>%
  left_join(
    baileydf %>%
      select(original_sample_id, Immune_Status)
  ) %>%
  left_join(
    mahapdf %>%
      select(original_sample_id, Immune_Status),
    by = "original_sample_id"
  ) %>%
  left_join(lukowskidf, by = c("sample_accession", "original_sample_id")) %>%
  mutate(temp = case_when(
    !is.na(Immune_Status.x) ~ Immune_Status.x,
    !is.na(Immune_Status.y) ~ Immune_Status.y,
    !is.na(immune_status) ~ immune_status,
    TRUE ~ NA
  )) %>%
  mutate(Immune_Status = case_when(
    temp == "CTR-IS" ~ "IS",
    temp == "IS-GVHD" ~ "IS",
    temp == "RTR-IS" ~ "IS",
    temp == "IC" ~ "IC",
    temp == "IS" ~ "IS", # this was missing before
    temp == "immunosuppressed" ~ "IS",
    temp == "immunocompetent" ~ "IC",
    study_name == "Chitsazzadeh_2016" ~ "IC",
    study_name == "Veenstra_2023" ~ "IC",
    TRUE ~ NA
  )) %>%
  rename(Immune_Status_Detail = temp) %>%
  select(-Immune_Status.x, -Immune_Status.y) %>%
  mutate(Immune_Status = factor(Immune_Status, levels = c("IC", "IS")))
write_csv(immunodf, "data/Immune_Status_CIBERSORTx.csv")

# immunodf <- immunodf %>% filter(study_name == "Bailey_2023")

immunodf %>% count(Immune_Status)
immunodf %>% count(Immune_Status_Detail)
immunodf %>% count(study_name, Immune_Status)

infoCols <- c("Mixture", "P-value", "Correlation", "RMSE", "Absolute score (sig.score)")
immuneCellTypeCols <- setdiff(colnames(immunedf), infoCols)

pvals <- list()
for (i in 1:length(immuneCellTypeCols)) {
  ct <- immuneCellTypeCols[i]
  tres <- wilcox.test(formula(paste0("`", ct, "`", "~ Immune_Status")), data=immunodf)
  print(paste(ct, ":", tres$p.value))
  pvals[[ct]] <- tres$p.value
}
pvals <- unlist(pvals)
pvalsAdj <- p.adjust(pvals, method = "holm")

pvals[c("B cells memory", "T cells CD8")]
pvalsAdj[c("B cells memory", "T cells CD8")]

bcellPlot <- immunodf %>%
  filter(!is.na(Immune_Status)) %>%
  ggplot(aes(Immune_Status, `B cells memory`)) +
  geom_boxplot(aes(fill = Immune_Status)) +
  guides(fill = "none") +
  labs(x = "", y = "Memory B Cells") +
  scale_fill_brewer(palette = "Pastel2") +
  theme(
    text = element_text(size = 14)
  )
wilcox.test(`B cells memory` ~ Immune_Status, data = immunodf)
cd8Plot <- immunodf %>%
  filter(!is.na(Immune_Status)) %>%
  ggplot(aes(Immune_Status, `T cells CD8`)) +
  geom_boxplot(aes(fill = Immune_Status)) +
  guides(fill = "none") +
  labs(x = "", y = "CD8 T Cells") +
  scale_fill_brewer(palette = "Pastel2") +
  theme(
    text = element_text(size = 14)
  )
wilcox.test(`T cells CD8` ~ Immune_Status, data = immunodf)

immunePlots <- bcellPlot | cd8Plot



# (tumorSubpopCompPlots | rdebCompPlots) / (fibSubpopCompPlots | immunePlots)

mainFigure <- (diffCondPlot | tskCondPlot | diffRdebPlot | tskRdebPlot) /
  (inflamFibCondPlot | mesenchFibCondPlot | bcellPlot | cd8Plot)

# immunodf %>%
#   ggplot(aes(Immune_Status, `Plasma cells`)) +
#   geom_boxplot()
# immunodf %>%
#   ggplot(aes(Immune_Status, `B cells naive`)) +
#   geom_boxplot()
# immunodf %>%
#   ggplot(aes(Immune_Status, `T cells follicular helper`)) +
#   geom_boxplot()
# immunodf %>%
#   ggplot(aes(Immune_Status, `Eosinophils`)) +
#   geom_boxplot()

## Save figures
ggsave(
  filename = file.path(figureDir, "Figure4_Boxplots.svg"),
  plot = mainFigure,
  width = 8,
  height = 8
)


ggsave(
  filename = file.path(figureDir, "Figure4_Tumor_Subpop_SCC_KA_Boxplots.svg"),
  plot = tumorSubpopCompPlots,
  width = 4,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_Fibroblast_Subpop_SCC_KA_Boxplots.svg"),
  plot = fibSubpopCompPlots,
  width = 4,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_TSK_InflamFib_Scatter.svg"),
  plot = tskInflamPlot,
  width = 4,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_Tumor_Subpop_RDEB_Boxplots.svg"),
  plot = rdebCompPlots,
  width = 4,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_Fibroblast_RDEB_Boxplot.svg"),
  plot = mesenchRdebPlot,
  width = 3,
  height = 4
)

ggsave(
  filename = file.path(figureDir, "Figure4_ImmuneCells_IC_IS_Boxplots.svg"),
  plot = immunePlots,
  width = 4,
  height = 4
)

