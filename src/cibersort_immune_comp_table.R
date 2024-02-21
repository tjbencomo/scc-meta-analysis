library(readr)
library(dplyr)
library(stringr)

# dataDir <- file.path("~", "Downloads", "scc-meta-analysis-data")
dataDir <- file.path("data")
lm22df <- read_csv(file.path(dataDir, "cibersort", "Limma_LM22_SCC_KA_Deconv.csv"))
# metdf <- read_csv(file.path(dataDir, "metadata_final_cohort.csv"))

immunodf <- read_csv(file.path(dataDir, "Immune_Status_CIBERSORTx.csv"))


# mahapdf <- read_csv(file.path(dataDir, "Mahapatra_Immune_Info.csv"))
# baileydf <- read_csv(file.path(dataDir, "Bailey_Supplemental_File1.csv"))
# baileydf$original_sample_id <- str_c("PRNJA844527", baileydf$SEQ_ID, sep="_")
# lukowskidf <- read_csv("data/E-MTAB-6430.csv")
# lukowskidf <- lukowskidf %>% 
#   rename(original_sample_id = sample_id, sample_accession = ena_accession)
# 
# 
# df <- lm22df %>%
#   inner_join(metdf, by = c("Mixture" = "sample_id"))
# 
# immunodf <- df %>%
#   left_join(
#     baileydf %>%
#       select(original_sample_id, Immune_Status)
#   ) %>%
#   left_join(
#     mahapdf %>%
#       select(original_sample_id, Immune_Status),
#     by = "original_sample_id"
#   ) %>%
#   left_join(lukowskidf, by = c("sample_accession", "original_sample_id")) %>%
#   mutate(temp = case_when(
#     !is.na(Immune_Status.x) ~ Immune_Status.x,
#     !is.na(Immune_Status.y) ~ Immune_Status.y,
#     !is.na(immune_status) ~ immune_status,
#     TRUE ~ NA
#   )) %>%
#   mutate(Immune_Status = case_when(
#     temp == "CTR-IS" ~ "IS",
#     temp == "IS-GVHD" ~ "IS",
#     temp == "RTR-IS" ~ "IS",
#     temp == "IC" ~ "IC",
#     temp == "IS" ~ "IS", # this was missing before
#     temp == "immunosuppressed" ~ "IS",
#     temp == "immunocompetent" ~ "IC",
#     study_name == "Chitsazzadeh_2016" ~ "IC",
#     study_name == "Veenstra_2023" ~ "IC",
#     TRUE ~ NA
#   )) %>%
#   rename(Immune_Status_Detail = temp) %>%
#   select(-Immune_Status.x, -Immune_Status.y) %>%
#   mutate(Immune_Status = factor(Immune_Status, levels = c("IC", "IS")))

# immunodf <- immunodf %>% filter(study_name == "Bailey_2023")

immunodf %>% count(Immune_Status)
immunodf %>% count(study_name, Immune_Status)

infoCols <- c("Mixture", "P-value", "Correlation", "RMSE", "Absolute score (sig.score)")
immuneCellTypeCols <- setdiff(colnames(lm22df), infoCols)

pvalList <- list()
icValList <- list()
isValList <- list()
for (i in 1:length(immuneCellTypeCols)) {
  ct <- immuneCellTypeCols[i]
  tres <- wilcox.test(formula(paste0("`", ct, "`", "~ Immune_Status")), data=immunodf)
  print(paste(ct, ":", tres$p.value))
  pvalList[[ct]] <- tres$p.value
  tmpdf <- immunodf %>% filter(Immune_Status == "IC")
  icValList[[ct]] <- median(tmpdf[[ct]])
  tmpdf <- immunodf %>% filter(Immune_Status == "IS")
  isValList[[ct]] <- median(tmpdf[[ct]])
}

testdf <- data.frame(
  celltype = names(pvalList),
  IC_median = unlist(icValList),
  IS_median = unlist(isValList),
  pval = unlist(pvalList)
) %>% as_tibble() %>% arrange(pval)
write_csv(testdf, file.path(dataDir, "cibersort", "LM22_IC_vs_IS_Comparison.csv"))
