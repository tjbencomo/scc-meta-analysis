library(readr)
library(dplyr)
library(stringr)


dataDir <- file.path("data")
decondf <- read_csv(file.path(dataDir, "cibersort", "Limma_NonImmune_SCC_KA_Deconv.csv"))
metdf <- read_csv(file.path(dataDir, "metadata_final_cohort.csv"))

df <- decondf %>%
  inner_join(metdf, by = c("Mixture" = "sample_id"))


df <- df %>%
  mutate(RDEB_Status = case_when(
    condition == "SCC" & study_name == "Cho_2018" ~ "RDEB",
    condition == "SCC" & study_name != "Cho_2018" ~ "Sporadic",
    TRUE ~ "Other"
  ))


tumorCellTypes <- c("Tumor_KC_Basal", "Tumor_KC_Cyc", "Tumor_KC_Diff", "TSK")
fibCellTypes <- c("Mesenchymal_Fibroblast", "Inflam_Fibroblast")
cellTypes <- c(tumorCellTypes, fibCellTypes)

rdebValList <- list()
sporadicValList <- list()
rdebPvalList <- list()
for (i in 1:length(cellTypes)) {
  ct <- cellTypes[i]
  tres <- wilcox.test(formula(paste0("`", ct, "`", "~ RDEB_Status")), data=df %>% filter(RDEB_Status != "Other"))
  print(paste(ct, ":", tres$p.value))
  rdebPvalList[[ct]] <- tres$p.value
  tmpdf <- df %>% filter(RDEB_Status == "RDEB")
  rdebValList[[ct]] <- median(tmpdf[[ct]])
  tmpdf <- df %>% filter(RDEB_Status == "Sporadic")
  sporadicValList[[ct]] <- median(tmpdf[[ct]])
}

condPvalList <- list()
sccValList <- list()
kaValList <- list()
for (i in 1:length(cellTypes)) {
  ct <- cellTypes[i]
  tres <- wilcox.test(formula(paste0("`", ct, "`", "~ condition")), data=df)
  print(paste(ct, ":", tres$p.value))
  condPvalList[[ct]] <- tres$p.value
  tmpdf <- df %>% filter(condition == "SCC")
  sccValList[[ct]] <- median(tmpdf[[ct]])
  tmpdf <- df %>% filter(condition == "KA")
  kaValList[[ct]] <- median(tmpdf[[ct]])
}

rdeb_df <- data.frame(
  comparison = "RDEB vs Sporadic",
  celltype = names(rdebPvalList),
  RDEB_median = unlist(rdebValList),
  Sporadic_median = unlist(sporadicValList),
  pval = unlist(rdebPvalList)
) %>% as_tibble()

write_csv(rdeb_df, file.path(dataDir, "cibersort", "NonImmune_RDEB_vs_Sporadic_Comparison.csv"))

cond_df <- data.frame(
  comparison = "KA vs SCC",
  celltype = names(condPvalList),
  SCC_median = unlist(sccValList),
  KA_median = unlist(kaValList),
  pval = unlist(condPvalList)
) %>% as_tibble()

write_csv(cond_df, file.path(dataDir, "cibersort", "NonImmune_SCC_vs_KA_Comparison.csv"))


