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
* `gsva_scoring.R` - calculate sample-wise enrichment of Hallmark and Reactome pathways with GSVA
* `investigate_outlier_expression.R` and `investigate_outliers_pca_umap.R` - investigate outliers
* `individual_study_deg_analysis.R` - perform DE analysis for NS vs SCC on each study cohort

## Figure Scripts
### Figure 1 and S1
* `pre_post_batch_correction_plots.R` - make PCA and UMAP plots
* `dvp_plots.R` - make DvP score overlay on PCA plot
* `studies_qc_scatterplot.R` - make QC plots in S1

### Figure 2 and S2
* `study_variation_plots.R`

### Figure 3 and S3
