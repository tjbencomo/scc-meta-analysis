library(readr)
library(dplyr)
library(fgsea)
library(ggplot2)
library(annotables)


plotRes <- function(geneset, ranks, title) {
  p <- plotEnrichment(geneset, ranks) +
    labs(
      x = "Rank", y = "Enrichment Score",
      title = title
    )
  
  p$layers[[1]]$aes_params$colour <- "black"
  p$layers[[1]]$aes_params$size <- 1.5
  p <- p + scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 0, colour = "white"),
      axis.ticks.x = element_line(size = 0),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 14, hjust = .5)
    )
  return(p)
}

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

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")

ns_ak_df <- read_csv(file.path(deseqDir, "AK_vs_NS.csv"))
ak_scc_df <- read_csv(file.path(deseqDir, "SCC_vs_AK.csv"))
ns_scc_df <- read_csv(file.path(deseqDir, "SCC_vs_NS.csv"))

# hnsc_aracne_df <- read_tsv("~/code/ret-r01/data/hnsc_aracne.tsv")
# retTargets <- hnsc_aracne_df %>% 
#   filter(Regulator == 5979, MoA > 0) %>%
#   inner_join(
#     grch38 %>% select(entrez, symbol),
#     by = c("Target" = "entrez")
#   ) %>% pull(symbol)
# write(retTargets, "data/HNSC_RET_Regulon.txt")
retTargets <- scan(file.path(dataDir, "HNSC_RET_Regulon.txt"), what=character())

ns_ak_ranks <- get_ranks(ns_ak_df)
ak_scc_ranks <- get_ranks(ak_scc_df)
ns_scc_ranks <- get_ranks(ns_scc_df)
gs <- list('RET_Targets' = retTargets)

set.seed(823)
ns_ak_res <- fgsea(gs, ns_ak_ranks, eps = 0)
ak_scc_res <- fgsea(gs, ak_scc_ranks, eps = 0)
ns_scc_res <- fgsea(gs, ns_scc_ranks, eps = 0)

ns_scc_plot <- plotRes(gs$RET_Targets, ns_scc_ranks, title = "HNSC RET Targets")
ns_scc_plot

ggsave(
  filename = file.path(figureDir, "FigureS3_NS_SCC_HNSC_RET_Targets.svg"),
  plot = ns_scc_plot,
  width = 4, height = 4
)
