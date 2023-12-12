# ## Description: Create mixture expression input file for CIBERSORTx with data from KA/SCC samples
# ## Needs to be in HGNC gene IDs just like the scRNA-seq data

library(DESeq2)
library(dplyr)
library(readr)
library(stringr)

dataDir <- file.path("data", "processed", "deseq")
outdir <- file.path("data", "cibersort")
outfp <- file.path(outdir, "limma_tumor_mixture_matrix.tsv")
samplesfp <- file.path(dataDir, "limma_batch_normalized.rds")

samples <- readRDS(samplesfp)
samples <- samples[, samples$condition %in% c("SCC", "KA")]
rownames(samples) <- rowData(samples)$symbol

M <- assay(samples)
idx <- rowMeans(M) > 3.3
M <- M[idx, ]
M[M < 0] <- 0
colnames(M) <- samples$sample_id
sampledf <- data.frame(sampleName = colnames(samples), sampleID = colnames(M))

write.table(M, file = outfp, sep = "\t", row.names = T, col.names = NA)

