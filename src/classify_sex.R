## Impute sex for each sample since not every paper provides sex information
## Doing GSVA scoring on VST data for male gene signature does not work
## XIST (raw counts) kind of works but not fully
## Use k-means clustering on XIST and chrY gene expression

library(DESeq2)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(annotables)
library(ggplot2)
library(patchwork)
library(plotly)

metadata <- read_csv("data/metadata_final_cohort.csv")

gs <- list('MaleGenes' = grch38 %>% filter(chr == "Y") %>% pull(symbol))

dds <- readRDS("data/processed/deseq/deseq_obj.rds")
rownames(dds) <- rowData(dds)$symbol

sexdf <- read_excel("data/Gender_Labels.xlsx")

metadata <- metadata %>%
  left_join(sexdf) %>%
  filter(condition != "BCC")

metadata %>%
  count(gender)

## This says XIST is only expressed in females and using expression of XIST is sufficient
## to classify samples in bulk data
## https://rdrr.io/github/Oshlack/speckle/man/classifySex.html
## Note batch correction + VST makes it so all samples have non-zero XIST expression
## Probably easier to just use raw counts - study noise seems unlikely to confound this

df <- data.frame(
  Sample = colnames(dds),
  XIST_norm = counts(dds, normalized = T)["XIST", ],
  XIST = counts(dds, normalized = F)["XIST", ]
)

## Score male genes
foundGenes <- gs$MaleGenes[gs$MaleGenes %in% rownames(dds)]
df$MaleScore <- colSums(counts(dds, normalized=F)[foundGenes, ])
df$MaleScoreNorm <- colSums(counts(dds, normalized=T)[foundGenes, ])

df <- df %>%
  inner_join(
    metadata %>%
      select(Sample, gender)
  )

df %>%
  ggplot(aes(gender, XIST_norm)) +
  geom_boxplot(aes(fill = gender))

df %>%
  ggplot(aes(gender, XIST)) +
  geom_boxplot(aes(fill = gender))


df %>%
  ggplot(aes(gender, MaleScore)) +
  geom_boxplot(aes(fill = gender))


df %>%
  ggplot(aes(gender, MaleScoreNorm)) +
  geom_boxplot(aes(fill = gender))

p1 <- df %>%
  ggplot(aes(log10(XIST), log10(MaleScore))) +
  geom_point(aes(color = gender))

p2 <- df %>%
  ggplot(aes(log10(XIST_norm), log10(MaleScoreNorm), label=Sample)) +
  geom_point(aes(color = gender), size = 3, alpha = .7) +
  labs(x = "Log10 XIST mRNA", y = "Log10 chrY Score") +
  theme_bw()
p2

p1 + p2


ggplotly(p2)

# Use k-means to group into male and female clusters
df$log_XIST_norm <- log1p(df$XIST_norm)
df$log_MaleScoreNorm <- log1p(df$MaleScoreNorm)
kmres <- kmeans(df[, c("log_XIST_norm", "log_MaleScoreNorm")], centers = 2, nstart = 100)
table(kmres$cluster)

df$km_cluster <- factor(kmres$cluster)
df %>%
  ggplot(aes(log10(XIST_norm), log10(MaleScoreNorm), label=Sample)) +
  geom_point(aes(color = km_cluster), size = 3, alpha = .7) +
  labs(x = "Log10 XIST mRNA", y = "Log10 chrY Score") +
  theme_bw()

df$inferred_sex <- ifelse(df$km_cluster == 1, "Female", "Male")

df %>%
  ggplot(aes(log10(XIST_norm), log10(MaleScoreNorm), label=Sample)) +
  geom_point(aes(color = inferred_sex), size = 3, alpha = .7) +
  labs(x = "Log10 XIST mRNA", y = "Log10 chrY Score") +
  theme_bw()

write_csv(
  df %>%
    rename(sex_from_study = gender) %>%
    select(Sample, log_XIST_norm, log_MaleScoreNorm, sex_from_study, inferred_sex),
  "data/inferred_sex_labels.csv"
)
