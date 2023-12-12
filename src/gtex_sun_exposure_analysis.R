## Description: Do sun exposed vs sun protected DE analysis on GTEx data.
## Identify Sun Exposed vs Protected (EvP) signature
##


library(DESeq2)
library(IHW)
library(readr)
library(dplyr)
library(BiocParallel)
register(MulticoreParam(8))

data_dir <- "data"
ddsfp <- file.path(data_dir, "GTEx_DESeq_object.rds")

if (!file.exists(ddsfp)) {
  print("Building DESeq2 object and running differential expression tests")
  countsfp <- file.path(data_dir, "Skin_GTex_Counts.csv")
  metadatafp <- file.path(data_dir, "GTEx_Skin_Samples2.csv")
  
  mat <- read_csv(countsfp)
  emat <- as.matrix(mat[, 3:dim(mat)[2]])
  rownames(emat) <- mat$Description
  colnames(emat) <- colnames(mat)[3:dim(mat)[2]]
  emat <- rowsum(emat, group = rownames(emat))
  rm(mat); gc()
  
  info <- read_csv(metadatafp)
  info$condition <- ifelse(info$SMTSD == "Skin - Not Sun Exposed (Suprapubic)", "NoSun", "Sun")
  info$condition <- factor(info$condition, levels = c("NoSun", "Sun"))
  info$sex <- factor(info$SEX)
  info$age <- factor(info$AGE)
  # info <- info[, c(1, 4)]
  
  stopifnot(all(colnames(emat) == info$SAMPID))
  dds <- DESeqDataSetFromMatrix(countData = emat,
                                colData = info,
                                design = ~ age + sex + condition)
  rm(emat); gc()
  dds <- DESeq(dds, parallel = TRUE)
  saveRDS(dds, ddsfp)
} else {
  print("Loading previously saved DESeq2 object")
  dds <- readRDS(ddsfp)
}
res <- results(dds, contrast = c("condition", "Sun", "NoSun"), filterFun = ihw, parallel=TRUE)
resMap <- lfcShrink(dds, coef = "condition_Sun_vs_NoSun", type = "apeglm", parallel=TRUE)

res_df <- data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(up_in_exposed = log2FoldChange > 0)
res_map_df <- data.frame(resMap) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(up_in_exposed = log2FoldChange > 0)

write_csv(res_df, file.path(data_dir, "Sun_Exposed_vs_Protected_DE_MLE.csv"))
write_csv(res_map_df, file.path(data_dir, "Sun_Exposed_vs_Protected_DE_MAP.csv"))
