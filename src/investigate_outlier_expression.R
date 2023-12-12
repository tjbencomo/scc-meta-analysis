## Description: Investigate expression of SCC markers and signature levels in outlier samples
## in Problem_Samples google sheets. Compare PTHLH, MMP1, DvP score, and E2F targets signature
## to help decide whether to exclude or swap sample.

library(DESeq2)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
dataDir <- file.path("data")
gseaDir <- file.path(dataDir, "gsea-results")
deseqDir <- file.path(dataDir, "processed", "deseq")

met_df <- read_csv(file.path(dataDir, "metadata_post_star_qc_cohort.csv"))
bailey_info <- read_csv(file.path(dataDir, "Bailey_Supplemental_File1.csv"))
dds <- readRDS(file.path(deseqDir, "deseq_obj.rds"))
rownames(dds) <- rowData(dds)$symbol
countMat <- counts(dds, normalize=T)

sig_df <- read_csv(file.path(gseaDir, "Hallmark_Reactome_DvP_GSVA_Scores.csv"))

rna_df <- data.frame(
  sample_id = dds$sample_id,
  condition = dds$condition,
  study = dds$study,
  original_sample_id = dds$original_sample_id,
  PTHLH = countMat["PTHLH", ],
  MMP1 = countMat["MMP1", ]
) %>% as_tibble()

rna_df <- rna_df %>%
  inner_join(
    sig_df %>% 
      mutate(sample_id = str_split(sampleID, "-", simplify=T)[, 1]) %>%
      select(sample_id, DvP_Score, HALLMARK_E2F_TARGETS),
    by = "sample_id"
  )

## Li 2020
li_df <- rna_df %>%
  filter(study == "Li_2020") %>%
  mutate(patient_id = str_split(original_sample_id, "_T|_N", simplify = T)[, 2]) %>%
  select(study, sample_id, original_sample_id, patient_id, condition, everything())

li_df %>%
  group_by(condition) %>%
  mutate(
    PTHLH_z = scale(log2(PTHLH))[, 1],
    MMP1_z = scale(log2(MMP1))[, 1],
    DvP_z = scale(DvP_Score)[, 1],
    E2F_z = scale(HALLMARK_E2F_TARGETS)[, 1]
  ) %>%
  View()

## Bailey 2023
bailey_df <- rna_df %>%
  filter(study == "Bailey_2023") %>%
  mutate(SEQ_ID = str_split(original_sample_id, "_", simplify = T)[, 2]) %>%
  inner_join(bailey_info)


