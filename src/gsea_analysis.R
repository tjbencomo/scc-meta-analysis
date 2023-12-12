## Description: Run GSEA analysis on genesets from MSigDB
## Save GSEA results files and draw plots for figures

library(readr)
library(readxl)
library(dplyr)
library(fgsea)
library(ggplot2)
library(patchwork)



get_ranks <- function(deg_df) {
  rank_df <- deg_df %>%
    filter(!is.na(gene_symbol), is.finite(stat)) %>%
    group_by(gene_symbol) %>%
    summarize(t_stat = mean(stat)) %>%
    ungroup()
  ranks <- rank_df$t_stat
  names(ranks) <- rank_df$gene_symbol
  return(ranks)
}

gsea_analysis <- function(ranks, gs, fp, up_colname) {
  gsea_res <- fgsea(gs, ranks, eps = 0, nPermSimple=10000)
  gsea_res[[up_colname]] <- gsea_res$NES > 0
  write_csv(gsea_res %>% select(-leadingEdge), fp)
}

dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")
outDir <- file.path(dataDir, "gsea-results")

# MSigDB genesets
reactome <- gmtPathways(file.path(dataDir, "genesets", "c2.cp.reactome.v7.4.symbols.gmt"))
hallmark <- gmtPathways(file.path(dataDir, "genesets", "h.all.v7.4.symbols.gmt"))
biocarta <- gmtPathways(file.path(dataDir, "genesets", "c2.cp.biocarta.v7.4.symbols.gmt"))


ns_ak_df <- read_csv(file.path(deseqDir, "AK_vs_NS.csv"))
ak_scc_df <- read_csv(file.path(deseqDir, "SCC_vs_AK.csv"))
ns_scc_df <- read_csv(file.path(deseqDir, "SCC_vs_NS.csv"))

ns_ak_ranks <- get_ranks(ns_ak_df)
ak_scc_ranks <- get_ranks(ak_scc_df)
ns_scc_ranks <- get_ranks(ns_scc_df)


# genesets: hallmark, reactome, biocarta
gs_list <- list(
  'Hallmark' = hallmark,
  'Reactome' = reactome,
  'Biocarta' = biocarta
)
set.seed(342435)
for (i in 1:length(gs_list)) {
  geneset_name <- names(gs_list)[i]
  print(paste("Running GSEA for NS vs AK with", geneset_name, "pathways"))
  fp <- file.path(outDir, paste0("NS_vs_AK_", geneset_name, "_GSEA.csv"))
  gsea_analysis(ns_ak_ranks, gs_list[[i]], fp, "up_in_AK")
}

set.seed(1324)
for (i in 1:length(gs_list)) {
  geneset_name <- names(gs_list)[i]
  print(paste("Running GSEA for AK vs SCC with", geneset_name, "pathways"))
  fp <- file.path(outDir, paste0("AK_vs_SCC_", geneset_name, "_GSEA.csv"))
  gsea_analysis(ak_scc_ranks, gs_list[[i]], fp, "up_in_SCC")
}

set.seed(18923)
for (i in 1:length(gs_list)) {
  geneset_name <- names(gs_list)[i]
  print(paste("Running GSEA for NS vs SCC with", geneset_name, "pathways"))
  fp <- file.path(outDir, paste0("NS_vs_SCC_", geneset_name, "_GSEA.csv"))
  gsea_analysis(ns_scc_ranks, gs_list[[i]], fp, "up_in_SCC")
}




