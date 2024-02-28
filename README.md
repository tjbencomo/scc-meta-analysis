# scc-meta-analysis
This repo contains code for "Gene expression landscape of cutaneous squamous cell carcinoma progression" by Bencomo and Lee. 

## Data
### Raw Data
Raw data (FASTQ files) can easily be downloaded from the apropriate sites with the data download scripts:

* `download_sra.sh` - download all FASTQ files available on SRA. Pulls SRR IDs from `srr_ids.txt`
* `download_ena.sh` - download all FASTQ files available on BioStudies. Pulls ERR IDs from `ena_fastq_urls.txt`

`srr_ids.txt` and `ena_fastq_urls.txt` can be found on Zenodo. Once the FASTQ files are downloaded, the Snakemake pipeline mentioned
above can generate BAM files and perform the DESeq2 analysis

### Processed Data
Raw FASTQ data were processed using a Snakemake workflow with STAR and RSEM with the hg38 (GENCODE v38) reference. 
It can be found [here](https://github.com/tjbencomo/nmsc-star).

**The processed data for this project (metadata, expression data, DEGs, fusions etc) is available at [Zenodo](https://zenodo.org/records/10272679).**

### `data/`
This folder contains various files needed to run some of the analysis scripts:

* `genesets/` - MSigDB geneset `.gmt` files
* `gsea-results/` - CSV file with results of GSVA scores for Hallmark and Reactome signatures on all samples
* `ji-2019/` - CSV file with fibroblast cell type subclassifications used for Ji data
* `qc/` - STAR and RSEM quality metrics for eac sample
* `schutz-2023` - Excel files from Schutz 2023 study with signatures for fibroblast subpopulations
* `Bailey_2023_Signatures.xlsx` and `Bailey_Supplemental_File1.csv` - metadata files from Bailey 2023 study
* `Up_Gene_Tally.csv` and `Dn_Gene_Tally.csv` - CSV files with number info on how many times each gene is up/downregulated in individually analyzed RNA-Seq studies (Figure 2)
* `Gender_Labels.xlsx` - sex/gender labels for samples reported in original studies. Some studies used the term sex while others used gender
* `Mahapatra_Immune_Info.csv` - Immune status info from Mahapatra 2020
* `inferred_sex_labels.csv` - sex information for all samples. Includes labels from original studies (`sex_from_study`), labels inferred using XIST and chrY expression (`inferred_sex`), and a combo of these two (`final_sex_label`) that uses `sex_from_study` if present, and otherwise uses `inferred_sex`
* `metadata_all_studies.csv` - sample info for all samples initially considered
* `metadata_post_star_qc_cohort.csv` - samples after filtering for STAR quality filters
* `metadata_final_cohort.csv` - final cohort of samples used in analyses


## Analysis Scripts
These scripts must be run before trying to generate figures

* `batch_correction.R` - perform limma batch correction on VST data
* `gsva_scoring.R` - calculate sample-wise enrichment of Hallmark and Reactome pathways with `GSVA`
* `investigate_outlier_expression.R` and `investigate_outliers_pca_umap.R` - investigate outliers
* `individual_study_deg_analysis.R` - perform DE analysis for NS vs SCC on each study cohort
* `gsea_analysis.R` - GSEA using Hallmark and Reactome pathways with `fgea`
* `Create_Ji_scRNA_SeuratObject.R` - Create Seurat object with Ji NS and SCC data
* `cibersort_mixture_maker.R` and `cibersort_scRNA_matrix_maker.R` - make CIBERSORTx input data
* `gtex_sun_exposure_analysis.R` - DE analysis of sun exposed vs protected skin from GTEx

### Fusion Calling
Fusions were detected with the STAR-Fusion pipeline. The workflow to run STAR-Fusion on the patient samples is in `workflow/fusions`.
`call_fusions_in_cell_lines.sh` runs STAR-Fusion on the cell line samples. 
Use `gather_fusions.py` to collate fusion results after STAR-Fusion finishes. 

## Figure Scripts
### Figure 1 and S1
* `pre_post_batch_correction_plots.R` - make PCA and UMAP plots
* `dvp_plots.R` - make DvP score overlay on PCA plot
* `studies_qc_scatterplot.R` - make QC plots in S1

### Figure 2 and S2
* `study_variation_plots.R`

### Figure 3, S3, and S4
* `volcano_plots.R` - draw DEG volcano plots in 3B and S3A
* `canonical_cscc_degs.R` - draw 3C heatmap showing canonical gene expression
* `ns_ak_scc_gsea_plots.R` - draw GSEA volcano plot in 3D
* `top_degs_in_scRNA.R` - single cell expression plots in 3D and S3E
* `highlighted_de_genes.R` - plot expression of RET ligand and targets in S3B
* `gsea_RET_signature.R` - draw GSEA plot of RET regulon
* `meta_driver_gene_expression.R` - draw expression heatmap of potential driver genes
* `iec_ka_ak_comparison.R` - Compare KA vs SCC GSVA scores and draw S4 boxplots

### Figure 4 and S5
* `fibroblast_subpop_analysis.R` - Define fibroblast subpopulations in Ji's data and draw S5A-B
* `cibersort_dvp_plots.R` - compare cell type abundances in KA vs SCC, RDEB vs sporadic, and IC vs IS. Figure 4 and S5C
* `cibersort_immune_comp_table.R` and `cibersort_tumor_comp_tables.R` - make tables 4, 5, and 7

### Figure S6
* `fusion_analysis.R`

### Figure 5
* `gtex_figures.R` - also defines EvP signature
