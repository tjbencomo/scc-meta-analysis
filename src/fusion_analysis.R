## Description: Analyze fusions identified by STAR-Fusion
## Compare with list of COSMIC driver genes and TCGA fusions from Gao 2018 Cell Reports paper

library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(UpSetR)
library(ggplot2)
library(patchwork)

# Set up paths
figureDir <- file.path("figures", "manuscript")
# dataDir <- file.path("data")
dataDir <- file.path("~", "Downloads", "scc-meta-analysis-data")

# Load fusions from patient samples
patient_fusions <- read_csv(file.path(dataDir, "fusions.csv"))
patient_fusions <- patient_fusions %>% 
  filter(condition != "BCC", condition != "AK_IEC", condition != "AK_IEC_SCC")
colnames(patient_fusions)[1] <- "FusionName"
patient_fusions$LeftGeneSymbol <- str_split(patient_fusions$LeftGene, "\\^", simplify = T)[, 1]
patient_fusions$RightGeneSymbol <- str_split(patient_fusions$RightGene, "\\^", simplify = T)[, 1]

# Load fusions from cell lines
cell_line_df <- read_csv(file.path(dataDir, "cell_line_fusions.csv"))
colnames(cell_line_df)[1] <- "FusionName"
cell_line_df$study_name <- "Cell Lines"
cell_line_df$LeftGeneSymbol <- str_split(cell_line_df$LeftGene, "\\^", simplify = T)[, 1]
cell_line_df$RightGeneSymbol <- str_split(cell_line_df$RightGene, "\\^", simplify = T)[, 1]

# Load COSMIC gene info
cosmicdf <- read_csv(file.path(dataDir, "cosmic_gene_info.csv"))

# Load fusions from TCGA Gao paper
tcgadf <- read_excel(file.path(dataDir, "tcga_fusions_gao.xlsx"), sheet = "Final fusion call set", skip = 1)
tcgadf$LeftGene <- str_split(tcgadf$Fusion, "--", simplify = T)[, 1]
tcgadf$RightGene <- str_split(tcgadf$Fusion, "--", simplify = T)[, 2]

# Load sample information
met_df <- read_csv(file.path(dataDir, "metadata_final_cohort.csv"))


# Load colors for graphics
condition_colors <- readRDS(file.path(dataDir, "condition_color_palette.rds"))
condition_colors <- c(condition_colors, c("A431" = "pink", "SCCIC1" = "orange"))

## Number of unique fusions in patient samples and cell lines
print(paste("Number of unique fusions in patient samples:", length(unique(patient_fusions$FusionName))))
print(paste("Number of unique fusions in cell lines:", length(unique(cell_line_df$FusionName))))

# Number of samples with a fusion detected
print(paste("Number of samples with at least one fusion:", length(unique(patient_fusions$sample_id))))
print(paste("Number of SCC with at least one fusion:", 
            length(unique(patient_fusions[patient_fusions$condition == "SCC", ]$sample_id))))

# Tally samples and adjust for 4 that crashed
samplesTested <- met_df %>%
  count(condition, name = "total_samples")
samplesTested[samplesTested$condition == "SCC", "total_samples"] <- 146 #1 crashed
samplesTested[samplesTested$condition == "AK", "total_samples"] <- 44 #2 crashed
samplesTested[samplesTested$condition == "AK_IEC", "total_samples"] <- 7 #1 crashed

fusion_frac_df <- patient_fusions %>%
  distinct(sample_id, condition) %>%
  count(condition, name = "n_samples") %>%
  inner_join(samplesTested) %>%
  mutate(Fusion_Frac = n_samples / total_samples)

condition_order <- c("NS", "AK", "IEC", "SCC", "KA")
fusionFracPlot <- fusion_frac_df %>%
  mutate(
    Fusion_Percent = Fusion_Frac * 100,
    condition = factor(condition, levels = condition_order)
  ) %>%
  ggplot(aes(condition, Fusion_Percent)) +
  geom_col(aes(fill = condition)) +
  scale_fill_manual(values = condition_colors) +
  theme_classic() +
  guides(fill = "none") +
  labs(x = "", y = "% With Fusions") +
  theme(text = element_text(size = 14))
fusionFracPlot

## Join sample fusions with cell line fusions
fusions <- patient_fusions %>%
  inner_join(met_df) %>%
  bind_rows(cell_line_df)


## Check for EGFR-PPARGC1A and ADCK4-NUMBL
egashiraFusions <- c("EGFR--PPARGC1A", "ADCK4--NUMBL")
fusions %>%
  filter(FusionName %in% egashiraFusions) %>%
  select(FusionName, LeftBreakpoint, RightBreakpoint, study_name, sample_id)


## Plot number of fusions per study and condition
fusionCounts <- fusions %>%
  count(study_name, sample_id, condition, name="n_fusions")
allSamples <- unique(c(met_df$sample_id, cell_line_df$sample_id))
zeroSamples <- allSamples[!(allSamples %in% fusionCounts$sample_id)]
missing_df <- data.frame(
  sample_id = zeroSamples,
  n_fusions = 0
) %>% left_join(met_df %>% select(sample_id, study_name, condition))

condition_order <- c("NS", "AK", "IEC", "SCC", "KA", "A431", "SCCIC1")
fusionCounts <- fusionCounts %>%
  bind_rows(missing_df) %>%
  filter(condition != "AK_IEC", condition != "AK_IEC_SCC", condition != "BCC") %>%
  mutate(condition = factor(condition, levels = condition_order))

studyCountPlot <- fusionCounts %>%
  ggplot(aes(reorder(study_name, n_fusions), n_fusions)) +
  geom_boxplot(aes(fill = study_name)) +
  theme_classic() +
  labs(x = "", y = "Fusions") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


conditionCountPlot <- fusionCounts %>%
  mutate(condition = factor(condition, )) %>%
  ggplot(aes(condition, n_fusions)) +
  geom_boxplot(aes(fill = condition)) +
  theme_classic() +
  labs(x = "", y = "Fusions") +
  guides(fill = "none") +
  scale_fill_manual(values = condition_colors)

lm(n_fusions ~ condition + study_name, data = fusionCounts) %>% summary()
lm(n_fusions ~ condition + study_name + 0, data = fusionCounts) %>% summary()

count_plots <- studyCountPlot | conditionCountPlot + plot_layout(widths = c(3, 1))
count_plots

## Get fusions found in each sample type
sccFusions <- fusions %>% filter(condition == "SCC") %>% pull(FusionName)
nsFusions <-  fusions %>% filter(condition == "NS") %>% pull(FusionName)
akFusions <-  fusions %>% filter(condition == "AK") %>% pull(FusionName)
# akIecFusions <- fusions %>% filter(condition == "AK_IEC") %>% pull(FusionName)
# akIecSccFusions <- fusions %>% filter(condition == "AK_IEC_SCC") %>% pull(FusionName)
iecFusions <-  fusions %>% filter(condition == "IEC") %>% pull(FusionName)
kaFusions <- fusions %>% filter(condition == "KA") %>% pull(FusionName)
a431Fusions <- fusions %>% filter(condition == "A431") %>% pull(FusionName)
sccic1Fusions <- fusions %>% filter(condition == "SCCIC1") %>% pull(FusionName)

## See fusion overlap across conditions and cell lines
fusion_list <- list(
  'SCC' = sccFusions,
  'NS' = nsFusions,
  'AK' = akFusions,
  # 'AK_IEC' = akIecFusions,
  # 'AK_IEC_SCC' = akIecSccFusions,
  'IEC' = iecFusions,
  'KA' = kaFusions,
  'A431' = a431Fusions,
  'SCCIC1' = sccic1Fusions
)
fusion_list <- sapply(fusion_list, unique)
patientSampleTypes <- unique(patient_fusions$condition)
fusion_list_no_cell_lines <- fusion_list[c("NS", "AK", "IEC", "SCC", "KA")]
patientUpsetPlot <- upset(
  fromList(fusion_list_no_cell_lines),
  order.by = "freq",
  text.scale = 2,
  nsets = length(fusion_list_no_cell_lines)
)
patientUpsetPlot

cellLinesUpsetPlot <- upset(
  fromList(fusion_list[c("SCC", "A431", "SCCIC1")]),
  order.by = "freq",
  text.scale = 2
)
cellLinesUpsetPlot

## Of all the unique fusions found in patient samples, how many are present in normal skin
print(paste(sum(unique(patient_fusions$FusionName) %in% fusion_list$NS), "of", 
            length(unique(patient_fusions$FusionName)), 
            "fusions from patient samples are found in normal skin samples"))

## How many fusions are only in lesional tissue?
lesionTypes <- c("AK", "IEC", "SCC", "KA")
lesionalFusions <-
  setdiff(
    unlist(fusion_list[lesionTypes]),
    fusion_list$NS
  )
print(paste(length(lesionalFusions), "fusions found that are only in lesional tissue and not in normal skin"))


## Describe characteristics of the lesional fusions
## How many in TCGA samples
## How many in normal tissue from other data banks
## How many have COSMIC gene(s)

lesionalFusionFeatures <- fusions %>%
  filter(condition %in% lesionTypes, FusionName %in% lesionalFusions) %>%
  distinct(FusionName, .keep_all = TRUE) %>%
  mutate(
    GTEx_Flag = str_detect(annots, "GTEx_recurrent_StarF2019"),
    Babiceanu_Normal_Flag = str_detect(annots, "Babiceanu_Normal"),
    Neighbors_Flag = str_detect(annots, "NEIGHBORS"),
    TCGA_Star_Flag = str_detect(annots, "TCGA_StarF2019"),
  ) %>%
  mutate(Normal_Tissue_Flag = case_when(
    GTEx_Flag | Babiceanu_Normal_Flag  ~ "Found in Normal Tissue",
    !GTEx_Flag & !Babiceanu_Normal_Flag ~ "Not Found in Normal Tissue",
    TRUE ~ "Something went wrong"
  )) %>%
  mutate(
    TCGA_Flag = case_when(
      TCGA_Star_Flag | (FusionName %in% tcgadf$Fusion) ~ "Found in TCGA Tumor(s)",
      TRUE ~ "Not Found in TCGA"
    )
  ) %>%
  mutate(
    LeftGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
    RightGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
  ) %>%
  mutate(Includes_COSMIC_Partner = LeftGene_Is_COSMIC | RightGene_Is_COSMIC)

propNormal <- sum(lesionalFusionFeatures$Normal_Tissue_Flag == "Found in Normal Tissue")
propTCGA <- sum(lesionalFusionFeatures$TCGA_Flag == "Found in TCGA Tumor(s)")
propCosmic <- sum(lesionalFusionFeatures$Includes_COSMIC_Partner == TRUE)

plot_df <- data.frame(
  variable = c("Found in Normal\nTissue", "Found in TCGA\nTumor(s)", "Has COSMIC\nGene(s)"),
  value = c(propNormal, propTCGA, propCosmic) / length(lesionalFusions) * 100
)

fusionInfoPlot <- plot_df %>%
  ggplot(aes(reorder(variable, -value), value)) +
  geom_col(aes(fill = variable), width = .7) +
  guides(fill = "none") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  labs(x = "", y = "% Fusions From Lesional Tissue") +
  scale_y_continuous(breaks = seq(0, 15, by = 5), limits = c(0, 15))
fusionInfoPlot


nonNormalLesionalFusions <- lesionalFusionFeatures %>%
  filter(
    Normal_Tissue_Flag != "Found in Normal Tissue"
  )
tcgaLesionalFusions <- lesionalFusionFeatures %>%
  filter(
    TCGA_Flag == "Found in TCGA Tumor(s)"
  )


nonNormalLesionFusionalTable <- fusions %>%
  filter(condition %in% lesionTypes, FusionName %in% nonNormalLesionalFusions$FusionName) %>%
  mutate(
    GTEx_Flag = str_detect(annots, "GTEx_recurrent_StarF2019"),
    Babiceanu_Normal_Flag = str_detect(annots, "Babiceanu_Normal"),
    Neighbors_Flag = str_detect(annots, "NEIGHBORS"),
    TCGA_Star_Flag = str_detect(annots, "TCGA_StarF2019"),
  ) %>%
  mutate(Normal_Tissue_Flag = case_when(
    GTEx_Flag | Babiceanu_Normal_Flag  ~ TRUE,
    !GTEx_Flag & !Babiceanu_Normal_Flag ~ FALSE,
  )) %>%
  mutate(
    TCGA_Flag = case_when(
      TCGA_Star_Flag | (FusionName %in% tcgadf$Fusion) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  mutate(
    LeftGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
    RightGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
  ) %>%
  mutate(Includes_COSMIC_Partner = LeftGene_Is_COSMIC | RightGene_Is_COSMIC) %>%
  count(FusionName, Normal_Tissue_Flag, TCGA_Flag, Includes_COSMIC_Partner, 
        LeftBreakpoint, RightBreakpoint, condition) %>%
  tidyr::pivot_wider(names_from = "condition", values_from = "n") %>%
  replace(is.na(.), 0) %>%
  select(FusionName, LeftBreakpoint, RightBreakpoint, SCC, AK, IEC, KA,
         Normal_Tissue_Flag, TCGA_Flag, Includes_COSMIC_Partner) %>% 
  arrange(desc(SCC))
ffpmdf <- fusions %>%
  filter(condition %in% lesionTypes, FusionName %in% nonNormalLesionalFusions$FusionName) %>%
  group_by(FusionName) %>%
  summarize(mean_FFPM = mean(FFPM))
xids <- paste0(nonNormalLesionFusionalTable$FusionName, "_",
              nonNormalLesionFusionalTable$LeftBreakpoint, "_",
              nonNormalLesionFusionalTable$RightBreakpoint)
proteinInfo <- fusions %>%
  mutate(fid = paste0(FusionName, "_", LeftBreakpoint, "_", RightBreakpoint)) %>%
  distinct(fid, .keep_all = T) %>%
  filter(fid %in% xids) %>%
  select(fid, PROT_FUSION_TYPE, FUSION_CDS, FUSION_TRANSL)

nonNormalLesionFusionalTable$fid <- xids
nonNormalLesionFusionalTable <- nonNormalLesionFusionalTable %>%
  inner_join(proteinInfo) %>% 
  select(-fid) %>%
  left_join(ffpmdf) %>%
  select(
    FusionName, LeftBreakpoint, RightBreakpoint, SCC, AK, IEC, KA, mean_FFPM, 
    Normal_Tissue_Flag, TCGA_Flag, Includes_COSMIC_Partner,
    PROT_FUSION_TYPE, FUSION_CDS, FUSION_TRANSL
  )

# tcgaLesionalFusionTable <- fusions %>%
#   filter(condition %in% lesionTypes, FusionName %in% tcgaLesionalFusions$FusionName) %>%
#   mutate(
#     GTEx_Flag = str_detect(annots, "GTEx_recurrent_StarF2019"),
#     Babiceanu_Normal_Flag = str_detect(annots, "Babiceanu_Normal"),
#     Neighbors_Flag = str_detect(annots, "NEIGHBORS"),
#     TCGA_Star_Flag = str_detect(annots, "TCGA_StarF2019"),
#   ) %>%
#   mutate(Normal_Tissue_Flag = case_when(
#     GTEx_Flag | Babiceanu_Normal_Flag  ~ TRUE,
#     !GTEx_Flag & !Babiceanu_Normal_Flag ~ FALSE,
#   )) %>%
#   mutate(
#     TCGA_Flag = case_when(
#       TCGA_Star_Flag | (FusionName %in% tcgadf$Fusion) ~ TRUE,
#       TRUE ~ FALSE
#     )
#   ) %>%
#   mutate(
#     LeftGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
#     RightGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
#   ) %>%
#   mutate(Includes_COSMIC_Partner = LeftGene_Is_COSMIC | RightGene_Is_COSMIC) %>%
#   count(FusionName, Normal_Tissue_Flag, TCGA_Flag, Includes_COSMIC_Partner, 
#         LeftBreakpoint, RightBreakpoint, condition) %>%
#   tidyr::pivot_wider(names_from = "condition", values_from = "n") %>%
#   replace(is.na(.), 0) %>%
#   select(FusionName, LeftBreakpoint, RightBreakpoint, SCC, AK, IEC, 
#          Normal_Tissue_Flag, TCGA_Flag, Includes_COSMIC_Partner) %>% 
#   arrange(desc(SCC))


write_csv(nonNormalLesionFusionalTable, file.path(dataDir, "Lesional_Fusions_Table1.csv"))
# write_csv(tcgaLesionalFusionTable, file.path(dataDir, "Lesional_Fusions_Table2.csv"))




# recurrence_stats <- fusions %>%
#   mutate(Present_TCGA = str_detect(annots, "TCGA_StarF2019")) %>%
#   mutate(Present_GTEx = str_detect(annots, "GTEx_recurrent_StarF2019")) %>%
#   mutate(Present_Normal = str_detect(annots, "Babiceanu_Normal") | Present_GTEx) %>%
#   count(FusionName, Present_TCGA, Present_GTEx, Present_Normal, LeftGeneSymbol, RightGeneSymbol, condition) %>%
#   tidyr::pivot_wider(names_from = "condition", values_from = "n") %>%
#   replace(is.na(.), 0) %>%
#   mutate(cosmic_partner = LeftGeneSymbol %in% cosmicdf$`Gene Symbol` | RightGeneSymbol %in% cosmicdf$`Gene Symbol`)
# 
# ## SCC specific fusions
# cancer_fusions <- recurrence_stats %>%
#   filter(NS == 0) %>%
#   arrange(desc(SCC))
# 
# ## KA specific fusions
# cancer_fusions %>%
#   filter(IEC > 0)
# 
# ## Save fusion table results
# cancer_fusions %>%
#   rename(Cosmic_Partner = cosmic_partner) %>%
#   select(-LeftGeneSymbol, -RightGeneSymbol) %>%
#   select(FusionName, Present_TCGA, Present_GTEx, Present_Normal, Cosmic_Partner, everything()) %>%
#   write_csv(file.path(dataDir, "lesional_fusions.csv"))
# 
# 
# ## Check Cosmic genes
# fusions %>%
#   filter(LeftGeneSymbol %in% cosmicdf$`Gene Symbol` | RightGeneSymbol %in% cosmicdf$`Gene Symbol`) %>%
#   count(condition)
#   
# goodFusions <- fusions %>%
#   filter(LeftGeneSymbol %in% cosmicdf$`Gene Symbol` | RightGeneSymbol %in% cosmicdf$`Gene Symbol`) %>%
#   filter(condition == "SCC", !(FusionName %in% c(nsFusions, akFusions, iecFusions)))
# 
# 
# ## Check TCGA Fusions
# sum(fusions$FusionName %in% tcgadf$Fusion)
# fusions %>%
#   filter(FusionName %in% tcgadf$Fusion) %>%
#   count(condition)
# 
# fusions %>%
#   filter(FusionName %in% tcgadf$Fusion, condition == "SCC")
# 
# 
# ## Check for Normal DB/TCGA/Cosmic genes/NEIGHBOR flags
# 
# ## Remove fusions found in normal samples from consideration
# ## Still consider fusions that are found in lesional samples that also appear
# ## in normal samples
# plot_info <- fusions %>%
#   filter(condition != "NS") %>%
#   distinct(FusionName, .keep_all = T) %>%
#   mutate(
#     GTEx_Flag = str_detect(annots, "GTEx_recurrent_StarF2019"),
#     Babiceanu_Normal_Flag = str_detect(annots, "Babiceanu_Normal"),
#     Neighbors_Flag = str_detect(annots, "NEIGHBORS"),
#     TCGA_Star_Flag = str_detect(annots, "TCGA_StarF2019"),
#     NS_Flag = FusionName %in% nsFusions
#   ) %>%
#   mutate(Normal_Tissue_Flag = case_when(
#     GTEx_Flag | Babiceanu_Normal_Flag | NS_Flag  ~ "Found in Normal Tissue",
#     !GTEx_Flag & !Babiceanu_Normal_Flag  & !NS_Flag ~ "Not Found in Normal Tissue",
#     TRUE ~ "Something went wrong"
#   )) %>%
#   mutate(
#     TCGA_Flag = case_when(
#       TCGA_Star_Flag | (FusionName %in% tcgadf$Fusion) ~ "Found in TCGA Tumor(s)",
#       TRUE ~ "Not Found in TCGA"
#     )
#   ) %>%
#   mutate(
#     LeftGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
#     RightGene_Is_COSMIC = LeftGeneSymbol %in% cosmicdf$`Gene Symbol`,
#   ) %>%
#   mutate(Includes_COSMIC_Partner = LeftGene_Is_COSMIC | RightGene_Is_COSMIC)
# 
# 
# propNormal <- sum(plot_info$Normal_Tissue_Flag == "Found in Normal Tissue")
# propTCGA <- sum(plot_info$TCGA_Flag == "Found in TCGA Tumor(s)")
# propCosmic <- sum(plot_info$Includes_COSMIC_Partner == TRUE)
# 
# plot_df <- data.frame(
#   variable = c("Found in Normal Tissue", "Found in TCGA Tumor(s)", "Has COSMIC Gene(s)"),
#   value = c(propNormal, propTCGA, propCosmic) / nrow(plot_info) * 100
# )
# 
# fusionInfoPlot <- plot_df %>%
#   ggplot(aes(reorder(variable, -value), value)) +
#   geom_col(aes(fill = variable), width = .7) +
#   guides(fill = "none") +
#   theme_classic() +
#   labs(x = "", y = "% Fusions From Lesional Tissue") +
#   theme(axis.text.x = element_text(angle = 35, hjust=1, vjust=1)) +
#   scale_y_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70))
# fusionInfoPlot


## Save plots
svg(file.path(figureDir, "Figure_S5_Conditions_Fusion_Overlap.svg"), width = 14, height = 8)
patientUpsetPlot
dev.off()

svg(file.path(figureDir, "Figure_S5_SCC_Cell_Lines_Fusion_Overlap.svg"), width = 10, height = 8)
cellLinesUpsetPlot
dev.off()

ggsave(
  file.path(figureDir, "FigureS5_Fusioned_Samples.svg"),
  plot = fusionFracPlot,
  width = 5, 
  height = 5
)

ggsave(
  filename = file.path(figureDir, "FigureS5_Fusion_Counts.svg"),
  plot = count_plots,
  width = 14,
  height = 6
)

ggsave(
  filename = file.path(figureDir, "FigureS5_Fusion_Attributes.svg"),
  plot = fusionInfoPlot,
  width = 5,
  height = 5
)
