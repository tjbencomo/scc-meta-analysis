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
above will generate BAM files and DESeq2 analysis

## Analysis Scripts
* `batch_correction.R` - perform limma batch correction on VST data
* `gsva_scoring.R` - calculate sample-wise enrichment of Hallmark and Reactome pathways with `GSVA`
* `investigate_outlier_expression.R` and `investigate_outliers_pca_umap.R` - investigate outliers
* `individual_study_deg_analysis.R` - perform DE analysis for NS vs SCC on each study cohort
* `gsea_analysis.R` - GSEA using Hallmark and Reactome pathways with `fgea`
* `Create_Ji_scRNA_SeuratObject.R` - Create Seurat object with Ji NS and SCC data

## Figure Scripts
### Figure 1 and S1
* `pre_post_batch_correction_plots.R` - make PCA and UMAP plots
* `dvp_plots.R` - make DvP score overlay on PCA plot
* `studies_qc_scatterplot.R` - make QC plots in S1

### Figure 2 and S2
* `study_variation_plots.R`

### Figure 3 and S3
* `volcano_plots.R` - draw DEG volcano plots in 3B and S3A
* `canonical_cscc_degs.R` - draw 3C heatmap showing canonical gene expression
* `ns_ak_scc_gsea_plots.R` - draw GSEA volcano plot in 3D
* `top_degs_in_scRNA.R` - single cell expression plots in 3D and S3E
* `highlighted_de_genes.R` - plot expression of RET ligand and targets in S3B
* `gsea_RET_signature.R` - draw GSEA plot of RET regulon
* `meta_driver_gene_expression.R` - draw expression heatmap of potential driver genes
