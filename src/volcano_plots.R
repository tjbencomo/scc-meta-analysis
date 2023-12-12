## Description: Create volcano and MA plots for contrasts:
## 1) SCC vs NS
## 2) SCC vs AK
## 3) AK vs NS
## 4) IEC vs NS
## 5) SCC vs IEC
## 6) IEC vs AK


library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

theme_set(theme_classic())

plotVolcano <- function(df, plotTitle = "", lfcThreshold = 1, padjThreshold = .05, lfcMax = 16,
                        genes = NULL, ngenes = 5, color_scheme) {
  if (is.null(genes)) {
    print("Selecting top genes")
    topGenes <- df %>%
      filter(!is.na(log2FoldChange), abs(log2FoldChange) < lfcMax) %>%
      group_by(log2FoldChange > 0) %>%
      slice_min(padj, n = ngenes) %>%
      ungroup() %>%
      pull(gene_symbol)
  } else {
    topGenes <- genes
  }
  print(topGenes)
  p <- df %>%
    mutate(colorLabel = case_when(
      log2FoldChange > lfcThreshold & padj < padjThreshold ~ "Up",
      log2FoldChange < -lfcThreshold & padj < padjThreshold ~ "Dn",
      TRUE ~ "Unchanged"
    )) %>%
    mutate(txtLabel = case_when(
      gene_symbol %in% topGenes ~ gene_symbol,
      TRUE ~ ""
    )) %>%
    filter(abs(log2FoldChange) < lfcMax) %>%
    # slice_sample(n = 1e3) %>%
    ggplot(aes(log2FoldChange, -log10(padj), label = txtLabel)) +
    geom_point(aes(color = colorLabel), alpha = .6) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", title = plotTitle) +
    geom_hline(yintercept = -log10(padjThreshold), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-lfcThreshold, lfcThreshold), linetype = "dashed", color = "red") +
    geom_text_repel(max.overlaps = Inf, min.segment.length = 0) +
    # scale_color_manual(values = c("Up" = "red", "Unchanged" = "gray", "Dn" = "blue")) +
    scale_color_manual(values = color_scheme) +
    theme(plot.title = element_text(hjust = .5, face = "bold"),
          text = element_text(size = 14)) +
    guides(color = "none")
  return(p)
}



data_dir <- file.path("data", "processed", "deseq")
figureDir <- file.path("figures", "manuscript")
up_down_colors <- readRDS("data/up_down_direction_colors.rds")
names(up_down_colors) <- c("Up", "Dn")
up_down_colors <- c(up_down_colors, c("Unchanged" = 'gray'))

padjThresh <- .01
lfcThres <- 10

scc_ns <- read_csv(file.path(data_dir, "SCC_vs_NS.csv"))
scc_ak <- read_csv(file.path(data_dir, "SCC_vs_AK.csv"))
ak_ns <- read_csv(file.path(data_dir, "AK_vs_NS.csv"))

## volcano plots
sccNsFdrThresh <- .001 # larger sample size = more stringent
sccNsPlot <- plotVolcano(scc_ns, plotTitle = "NS vs SCC", padjThreshold = sccNsFdrThresh, lfcMax = lfcThres, color_scheme = up_down_colors)
sccAkPlot <- plotVolcano(scc_ak, plotTitle = "AK vs SCC", padjThreshold = padjThresh, lfcMax = lfcThres, color_scheme = up_down_colors)
akNsPlot <- plotVolcano(ak_ns, plotTitle = "NS vs AK", padjThreshold = padjThresh, lfcMax = lfcThres, color_scheme = up_down_colors)


mainPlot <- akNsPlot | sccAkPlot
mainPlot

## Save plots
ggsave(
  filename = file.path(figureDir, "Figure3_Volcano_Plots_Stages.png"),
  plot = akNsPlot | sccAkPlot,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(figureDir, "Figure3_Volcano_Plots_NS_SCC.png"),
  plot = sccNsPlot,
  width = 4,
  height = 4
)
