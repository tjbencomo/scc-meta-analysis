## Description: Remove study-specific batch effects from RNA-seq data
## to avoid batch effects in ARACNe network generation

library(DESeq2)

dataDir <- file.path("data", "processed", "deseq")
vsd <- readRDS(file.path(dataDir, "vst_normalized_counts.rds"))

plotPCA(vsd, intgroup = "condition")
plotPCA(vsd, intgroup = "study")
mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$study, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "study")
plotPCA(vsd, intgroup = "condition")

saveRDS(vsd, file.path(dataDir, "limma_batch_normalized.rds"))
