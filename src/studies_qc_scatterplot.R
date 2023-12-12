## Description: Plot percent mapped vs total number of mapped reads for all studies with
## data available (excluding Arron)


library(readr)
library(dplyr)
library(ggplot2)

figureDir <- file.path("figures", "manuscript")
dataDir <- "data"
starDir <- file.path(dataDir, "qc", "star-rsem", "multiqc_data")
stardf <- read_tsv(file.path(starDir, "multiqc_star.txt"))
rsemdf <- read_tsv(file.path(starDir, "multiqc_rsem.txt"))
metadata <- read_csv(file.path(dataDir, "metadata_all_studies.csv"))
metadata$Sample <- paste0(metadata$sample_id, "-", metadata$condition)
# qualimapdf$Sample <- basename(dirname(qualimapdf$Sample))


# blacklist_studies <- c("Wan_2019")
blacklist_studies <- c()
df <- stardf %>% inner_join(metadata) %>%
  inner_join(rsemdf) %>%
  mutate(
    total_mapped = uniquely_mapped + multimapped,
    total_mapped_percent = uniquely_mapped_percent + multimapped_percent
  ) %>%
  filter(!(study_name %in% blacklist_studies))

qcPlot <- df %>%
  mutate(study_name = str_replace(study_name, "_", " ")) %>%
  ggplot(aes(total_mapped_percent, total_mapped)) +
  geom_point(aes(fill = study_name), size=2, pch=21) +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept = 5e6, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 60, color = "red", linetype = "dashed") +
  labs(x = "Percent Mapped Reads", y = "Log10 [Number of Mapped Reads]",
       fill = "Study") +
  scale_fill_brewer(palette = "Set3") +
  theme(
    text = element_text(size = 14)
  )
qcPlot


readTotalPlot <- df %>%
  mutate(reads_million = total_mapped / 1e6) %>%
  mutate(study_name = str_replace(study_name, "_", " ")) %>%
  ggplot(aes(reorder(study_name, reads_million), reads_million)) +
  geom_boxplot(aes(fill = study_name)) +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = "none") +
  labs(x = "", y = "Total Mapped Reads (Millions)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    text = element_text(size = 14)
  )
readTotalPlot


# Samples removed due to abberant expression in outlier analysis
# For stats described in figure 1 results section
blacklist <- c(
  "Sample_290",
  "Sample_261",
  "Sample_336",
  "Sample_53",
  "Sample_54",
  "Sample_57"
)
df %>%
  filter(study_name != "Wan_2019") %>%
  filter(!(sample_id %in% blacklist)) %>%
  mutate(reads_million = total_mapped / 1e6) %>%
  summarize(
    mean_reads = mean(reads_million),
    min_reads = min(reads_million),
    max_reads = max(reads_million)
  )

ggsave(
  filename = file.path(figureDir, "Figure1_QC_Scatter.svg"),
  plot = qcPlot,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(figureDir, "Figure 1_Study_Read_Total.svg"),
  plot = readTotalPlot,
  width = 4,
  height = 5
)
