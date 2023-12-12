## Description: Create signature matrix from Andrew scRNA-seq data for CIBERSORTx
## SCC deconvolution
## DC, B Cell, Cd4/cd8/Treg, Endothelial, Inflam/Mesench fib, LC, Mac, 
## MDSC, Melanocyte, PDC, Tumor basal/cyc/diff/TSK

library(Seurat)
library(readr)
library(dplyr)
library(stringr)

scDataDir <- file.path("data", "ji-2019")
outdir <- file.path("data", "cibersort")
full_tumor_outfp <- file.path(outdir, "scRNA_tumor_signature_matrix.tsv")
full_tumor_merged_outfp <- file.path(outdir, "scRNA_tumor_merge_labels.tsv")
tumor_cell_outfp <- file.path(outdir, "scRNA_tumor_cell_signature_matrix.tsv")
non_immune_outfp <- file.path(outdir, "scRNA_nonImmune_signature_matrix.tsv")

cellsfp <- file.path(scDataDir, "Ji_2019_Cells.rds")
fibroblastfp <- file.path(scDataDir, "fibroblast_annotations.csv")

cells <- readRDS(cellsfp)
mig_dc_types <- c("Mig_CD1C", "Mig_CLEC9A", "Mig_LC")
keepCells <- cells@meta.data %>%
  filter(
    tum.norm == "Tumor",
    level2_celltype != "Keratinocyte",
    level3_celltype != "Eccrine",
    level3_celltype != "Pilosebaceous",
    # level3_celltype != "MDSC",
    level3_celltype != "NK",
    !(level3_celltype %in% mig_dc_types),
  ) %>%
  tibble::rownames_to_column("barcode") %>%
  pull(barcode)
cells <- cells[, keepCells]

fib_df <- read_csv(fibroblastfp)

## Add fibroblast labels
## Collapse other labels - CD4 and CD8 T cells?
cells@meta.data <- cells@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(fib_df %>% select(barcode, fibroblast_type)) %>%
  mutate(fibroblast_type = str_replace(fibroblast_type, " ", "_")) %>%
  tibble::column_to_rownames("barcode")

cd4_tcells <- c("CD4_Exh", "CD4_Naive", "CD4_Pre_Exh", "CD4_RGCC")
cd8_tcells <- c("CD8_EM", "CD8_EMRA", "CD8_Exh", "CD8_Naive")
dendritic_cells <- c("ASDC", "CD1C", "CLEC9A")
normal_kc <- c("Normal_KC_Basal", "Normal_KC_Cyc", "Normal_KC_Diff")
labels <- cells@meta.data %>%
  mutate(label = case_when(
    level3_celltype == "Endothelial Cell" ~ "Endothelial",
    level3_celltype == "B Cell" ~ "B_Cell",
    level3_celltype %in% dendritic_cells ~ "DC",
    level3_celltype == "Fibroblast" ~ fibroblast_type,
    level3_celltype %in% cd4_tcells ~ "CD4_TCell",
    level3_celltype %in% cd8_tcells ~ "CD8_TCell",
    level3_celltype %in% normal_kc ~ "Normal_KC",
    TRUE ~ level3_celltype
  )) %>%
  pull(label)
table(labels)
cells$label <- labels

# See note 5 in the book chapter about how many cells/samples are needed
# for accurate signature construction
sigMat <- cells@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  group_by(patient, label) %>%
  slice_head(n = 12) %>%
  group_by(label) %>%
  mutate(cellnum = row_number()) %>%
  mutate(cellid = str_c(label, cellnum, sep = "."))
table(sigMat$label)

M <- cells@assays$RNA@counts
zeroGenes <- rowSums(M)
keepGenes <- zeroGenes[zeroGenes >5] %>% names()
M <- M[keepGenes, sigMat$barcode]
colnames(M) <- sigMat$label
# M <- as.data.frame(M)
# M$gene <- rownames(M)
# M <- M %>% select(gene, everything())

write.table(M, file = full_tumor_outfp, sep = "\t", row.names = T, col.names = NA)
merged_df <- data.frame(
  first_label = colnames(M)
) %>%
  mutate(second_label = case_when(
    first_label %in% c("B_Cell", "CD4_TCell", "CD8_TCell", "Treg") ~ "Lymphocytes",
    first_label %in% c("DC", "LC", "Mac", "PDC", "MDSC") ~ "Myeloid",
    first_label %in% c("Inflam_Fibroblast", "Mesenchymal_Fibroblast") ~ "Fibroblast",
    TRUE ~ first_label
  )) %>%
  select(second_label)
cat(paste(merged_df$second_label, collapse = "\t"), 
    file = full_tumor_merged_outfp, append = FALSE, sep = "\t")

## Only use tumor cell subpops for signature matrix
tumor_celltypes <- c("Tumor_KC_Basal", "Tumor_KC_Cyc", "Tumor_KC_Diff", "TSK")
# fibroblast_celltypes <- c("Inflam_Fibroblast", "Mesenchymal_Fibroblast")
# nonImmuneCellTypes <- c(tumor_celltypes, fibroblast_celltypes)
tumorSigMat <- cells@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  filter(tum.norm == "Tumor", level3_celltype %in% tumor_celltypes) %>%
  group_by(patient, label) %>%
  slice_head(n = 20) %>%
  group_by(label) %>%
  mutate(cellnum = row_number()) %>%
  mutate(cellid = str_c(label, cellnum, sep = "."))
table(tumorSigMat$label)

tumors <- subset(cells, tum.norm == "Tumor")

M_tumor <- tumors@assays$RNA@counts
tumorZeroGenes <- rowSums(M_tumor)
tumorKeepGenes <- tumorZeroGenes[tumorZeroGenes > 5] %>% names()
M_tumor <- M_tumor[tumorKeepGenes, tumorSigMat$barcode]
colnames(M_tumor) <- tumorSigMat$label

write.table(M_tumor, file = tumor_cell_outfp, sep = "\t", row.names = T, col.names = NA)

## Only use non-immune cell types for signature matrix
fibroblast_celltypes <- c("Inflam_Fibroblast", "Mesenchymal_Fibroblast")
other_celltypes <- c("Endothelial", "Melanocyte")
nonImmuneCellTypes <- c(tumor_celltypes, fibroblast_celltypes, other_celltypes)
nonImmuneSigMat <- cells@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  filter(tum.norm == "Tumor", label %in% nonImmuneCellTypes) %>%
  group_by(patient, label) %>%
  slice_head(n = 20) %>%
  group_by(label) %>%
  mutate(cellnum = row_number()) %>%
  mutate(cellid = str_c(label, cellnum, sep = "."))
table(nonImmuneSigMat$label)


M_nonimmune <- tumors@assays$RNA@counts
nonImmuneZeroGenes <- rowSums(M_nonimmune)
nonImmuneKeepGenes <- nonImmuneZeroGenes[nonImmuneZeroGenes > 5] %>% names()
M_nonimmune <- M_nonimmune[nonImmuneKeepGenes, nonImmuneSigMat$barcode]
colnames(M_nonimmune) <- nonImmuneSigMat$label

write.table(M_nonimmune, file = non_immune_outfp, sep = "\t", row.names = T, col.names = NA)
