library(DESeq2)
library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(theme_classic())

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data")
deseqDir <- file.path(dataDir, "processed", "deseq")

vsd <- readRDS(file.path(deseqDir, "limma_normalized_vst.rds"))
genes <- c("ARTN", "DUSP6", "SPRY4")
idx <- which(rowData(vsd)$symbol %in% genes)

rnadf <- data.frame(
  sampleID = colnames(vsd),
  condition = vsd$condition,
  t(assay(vsd)[idx, ])
) %>% filter(condition != "AK_IEC", condition != "AK_IEC_SCC") %>% as_tibble()
colnames(rnadf)[3:ncol(rnadf)] <- rowData(vsd)$symbol[idx]
rnadf$condition <- factor(rnadf$condition, levels = c("NS", "AK", "IEC", "SCC", "KA"))
condition_colors <- readRDS("data/condition_color_palette.rds")



artnPlot <- rnadf %>%
  ggplot(aes(condition, ARTN)) +
  geom_boxplot(aes(fill = condition)) +
  labs(x = "", y = "Log2 ARTN mRNA", title = "ARTN") +
  theme(plot.title = element_text(hjust = .5)) +
  guides(fill = "none") +
  scale_fill_manual(values = condition_colors)
artnPlot



dusp6Plot <- rnadf %>%
  ggplot(aes(condition, DUSP6)) +
  geom_boxplot(aes(fill = condition)) +
  labs(x = "", y = "Log2 DUSP6 mRNA", title = "DUSP6") +
  theme(plot.title = element_text(hjust = .5)) +
  guides(fill = "none") +
  scale_fill_manual(values = condition_colors)
dusp6Plot



spry4Plot <- rnadf %>%
  ggplot(aes(condition, SPRY4)) +
  geom_boxplot(aes(fill = condition)) +
  labs(x = "", y = "Log2 SPRY4 mRNA", title = "SPRY4") +
  theme(plot.title = element_text(hjust = .5)) +
  guides(fill = "none") +
  scale_fill_manual(values = condition_colors)
spry4Plot

p <- artnPlot | dusp6Plot | spry4Plot
p

ggsave(
  filename = file.path(figureDir, "FigureS3_RET_Genes.svg"),
  plot = p,
  width = 7,
  height = 4
)


