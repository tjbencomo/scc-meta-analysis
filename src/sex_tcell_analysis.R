## Analyze impact of sex on CD8 T cell abundance
## Use 3 different sex classifications:
## 1) Sex labels provided by original studies
## 2) Sex labels inferred from RNA data (XIST and chrY expression)
## 3) Combination of 1) and 2) - use 1) if exists otherwise use 2)


library(readr)
library(dplyr)
library(ggplot2)

decondf <- read_csv("data/Immune_Status_CIBERSORTx.csv")
sexdf <- read_csv("data/inferred_sex_labels.csv")

df <- decondf %>%
  inner_join(sexdf, by = c("Sample")) %>%
  filter(!is.na(Immune_Status))

# Abundance of each sex per immune status
df %>%
  count(Immune_Status, sex_from_study)
df %>%
  ggplot(aes(sex_from_study)) +
  geom_bar(aes(fill = sex_from_study)) +
  facet_wrap(~Immune_Status)

df %>%
  count(Immune_Status, inferred_sex)
df %>%
  ggplot(aes(inferred_sex)) +
  geom_bar(aes(fill = inferred_sex)) +
  facet_wrap(~Immune_Status)


df %>%
  count(Immune_Status, final_sex_label)
df %>%
  ggplot(aes(final_sex_label)) +
  geom_bar(aes(fill = final_sex_label)) +
  facet_wrap(~Immune_Status)



## Note these are only for SCC/KA

## Sex labels from studies
## How many of each sex
df %>%
  count(sex_from_study)

## Test if difference in proportion on males vs females + stratify by Immune Status
sexFromstudyImmuneMat <- table(df$Immune_Status, df$sex_from_study)
sexFromstudyImmuneMat
fisher.test(sexFromstudyImmuneMat)

df %>%
  ggplot(aes(sex_from_study, `T cells CD8`)) +
  geom_boxplot(aes(fill = sex_from_study))
df %>%
  ggplot(aes(sex_from_study, `T cells CD8`)) +
  geom_boxplot(aes(fill = sex_from_study)) +
  facet_wrap(~Immune_Status)
wilcox.test(`T cells CD8` ~ sex_from_study, data=df %>% filter(!is.na(sex_from_study)))
wilcox.test(`T cells CD8` ~ sex_from_study, data=df %>% filter(!is.na(sex_from_study), Immune_Status == "IC"))
wilcox.test(`T cells CD8` ~ sex_from_study, data=df %>% filter(!is.na(sex_from_study), Immune_Status == "IS"))



## Inferred sex from XIST + chrY gene expression
## How many of each sex
df %>%
  count(inferred_sex)

## Test if difference in proportion on males vs females + stratify by Immune Status
sexFromRNAImmuneMat <- table(df$Immune_Status, df$inferred_sex)
sexFromRNAImmuneMat
fisher.test(sexFromRNAImmuneMat)

df %>%
  ggplot(aes(inferred_sex, `T cells CD8`)) +
  geom_boxplot(aes(fill = inferred_sex))
df %>%
  ggplot(aes(inferred_sex, `T cells CD8`)) +
  geom_boxplot(aes(fill = inferred_sex)) +
  facet_wrap(~Immune_Status)
wilcox.test(`T cells CD8` ~ inferred_sex, data=df %>% filter(!is.na(inferred_sex)))
wilcox.test(`T cells CD8` ~ inferred_sex, data=df %>% filter(!is.na(inferred_sex), Immune_Status == "IC"))
wilcox.test(`T cells CD8` ~ inferred_sex, data=df %>% filter(!is.na(inferred_sex), Immune_Status == "IS"))

## Combination of sex from study + inferred sex for unknown samples
## How many of each sex
df %>%
  count(final_sex_label)

## Test if difference in proportion on males vs females + stratify by Immune Status
sexFinalImmuneMat <- table(df$Immune_Status, df$final_sex_label)
sexFinalImmuneMat
fisher.test(sexFinalImmuneMat)

df %>%
  ggplot(aes(final_sex_label, `T cells CD8`)) +
  geom_boxplot(aes(fill = final_sex_label))
df %>%
  ggplot(aes(final_sex_label, `T cells CD8`)) +
  geom_boxplot(aes(fill = final_sex_label)) +
  facet_wrap(~Immune_Status)
wilcox.test(`T cells CD8` ~ final_sex_label, data=df %>% filter(!is.na(final_sex_label)))
wilcox.test(`T cells CD8` ~ final_sex_label, data=df %>% filter(!is.na(final_sex_label), Immune_Status == "IC"))
wilcox.test(`T cells CD8` ~ final_sex_label, data=df %>% filter(!is.na(final_sex_label), Immune_Status == "IS"))

wilcox.test(`T cells CD8` ~ Immune_Status, data=df)
