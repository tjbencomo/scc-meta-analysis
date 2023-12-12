## Description: Analyze GSVA scores for Hallmark/Reactome genesets
## to analyze molecular changes inIEC to KA/SCC

library(RegParallel)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

figureDir <- file.path("figures", "manuscript")
dataDir <- file.path("data", "gsea-results")

gsdf <- read_csv(file.path(dataDir, "Hallmark_Reactome_DvP_GSVA_Scores.csv"))
gsdf <- gsdf %>% filter(condition != "AK_IEC", condition != "AK_IEC_SCC")
gsdf$condition <- factor(gsdf$condition, levels = c("NS", "AK", "IEC", "SCC", "KA"))
pathways <- colnames(gsdf)
pathways <- pathways[str_detect(pathways, "HALLMARK|REACTOME")]
hallmark <- pathways[str_detect(pathways, "HALLMARK")]
reactome <- pathways[str_detect(pathways, "REACTOME")]



ak_iec <- RegParallel(
  data = gsdf %>% filter(condition %in% c("AK", "IEC")),
  formula = '[*] ~ condition',
  FUN = function(formula, data)
    lm(formula = formula,
        data = data),
  FUNtype = 'lm',
  variables = pathways,
  p.adjust = "fdr"
)

iec_ka <- RegParallel(
  data = gsdf %>% filter(condition %in% c("KA", "IEC")),
  formula = '[*] ~ condition',
  FUN = function(formula, data)
    lm(formula = formula,
       data = data),
  FUNtype = 'lm',
  variables = pathways,
  p.adjust = "fdr"
)

scc_iec <- RegParallel(
  data = gsdf %>% filter(condition %in% c("SCC", "IEC")),
  formula = '[*] ~ condition',
  FUN = function(formula, data)
    lm(formula = formula,
       data = data),
  FUNtype = 'lm',
  variables = pathways,
  p.adjust = "fdr"
)

scc_ka <- RegParallel(
  data = gsdf %>% filter(condition %in% c("SCC", "KA")),
  formula = '[*] ~ condition',
  FUN = function(formula, data)
    lm(formula = formula,
       data = data),
  FUNtype = 'lm',
  variables = pathways,
  p.adjust = "fdr"
)


## Save results
# write_csv(ak_iec, file.path(dataDir, "AK_vs_IEC_GSVA.csv.gz"))
write_csv(iec_ka, file.path(dataDir, "IEC_vs_KA_GSVA.csv.gz"))
write_csv(scc_iec, file.path(dataDir, "IEC_vs_SCC_GSVA.csv.gz"))
write_csv(scc_ka, file.path(dataDir, "KA_vs_SCC_GSVA.csv.gz"))


## Plots
theme_set(theme_classic())
condition_colors <- readRDS("data/condition_color_palette.rds")
fdr_thresh <- .01
iec_ka %>%
  filter(Variable %in% hallmark) %>%
  mutate(color_label = case_when(
    P.adjust < fdr_thresh & Beta > 0 ~ "Up",
    P.adjust < fdr_thresh & Beta < 0 ~ "Down",
    TRUE ~ "ns"
  )) %>%
  ggplot(aes(Beta, -log10(P.adjust))) +
    geom_point(aes(fill = color_label), pch=21, size=2) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "ns" = "grey")) +
  labs(x = "Log2FC", y = "")

scc_iec %>%
  filter(Variable %in% hallmark) %>%
  mutate(color_label = case_when(
    P.adjust < fdr_thresh & Beta > 0 ~ "Up",
    P.adjust < fdr_thresh & Beta < 0 ~ "Down",
    TRUE ~ "ns"
  )) %>%
  ggplot(aes(Beta, -log10(P.adjust))) +
  geom_point(aes(fill = color_label), pch=21, size=2) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "ns" = "grey"))


notchPlot <- gsdf %>%
  ggplot(aes(condition, HALLMARK_NOTCH_SIGNALING)) +
  geom_boxplot(aes(fill = condition)) +
  labs(title = "NOTCH Signaling", x = "", y = "Hallmark NOTCH Signaling") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = condition_colors)
notchPlot
lm(HALLMARK_NOTCH_SIGNALING ~ condition, data=gsdf) %>% summary()

e2fPlot <- gsdf %>%
  ggplot(aes(condition, HALLMARK_E2F_TARGETS)) +
  geom_boxplot(aes(fill = condition)) +
  labs(title = "E2F Targets", x = "", y = "Hallmark E2F Targets") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = condition_colors)
e2fPlot
lm(HALLMARK_E2F_TARGETS ~ condition, data=gsdf) %>% summary()

mycPlot <- gsdf %>%
  ggplot(aes(condition, HALLMARK_MYC_TARGETS_V2)) +
  geom_boxplot(aes(fill = condition)) +
  labs(title = "MYC V2 Targets", x = "", y = "Hallmark MYC V2 Targets") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_fill_manual(values = condition_colors)
mycPlot
lm(HALLMARK_MYC_TARGETS_V2 ~ condition, data=gsdf) %>% summary()

p1 <- (e2fPlot | notchPlot | mycPlot) + plot_layout(guides = "collect")
p1


## Save plots
ggsave(
  filename = file.path(figureDir, "FigureS3_IEC_KA_SCC_Pathways.svg"),
  plot = p1,
  width = 10,
  height = 4
)

