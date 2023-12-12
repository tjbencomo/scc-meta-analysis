# scc-meta-analysis
This repo contains code for "Gene expression landscape of cutaneous squamous cell carcinoma progression" by Bencomo and Lee. 

## Data
### Processed Data
Raw FASTQ data were processed using a Snakemake workflow with STAR and RSEM. 
It can be found [here](https://github.com/tjbencomo/nmsc-star).
The processed data for this project is available at [Zenodo](https://zenodo.org/records/10272679).

### Raw Data
Raw data can easily be downloaded from the apropriate sites with the data download scripts:

* `download_sra.sh` - download all FASTQ files available on SRA. Pulls SRR IDs from `srr_ids.txt`
* `download_ena.sh` - download all FASTQ files available on BioStudies. Pulls ERR IDs from `ena_fastq_urls.txt`

`srr_ids.txt` and `ena_fastq_urls.txt` can be found on Zenodo. Once the FASTQ files are downloaded, the Snakemake pipeline mentioned
above can generate BAM files and perform the DESeq2 analysis

## Analysis Scripts
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
