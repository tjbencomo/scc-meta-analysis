# scc-meta-analysis
This repo contains code for "Gene expression landscape of cutaneous squamous cell carcinoma progression" by Bencomo and Lee. 

## Data Processing
Raw FASTQ data were processed using a Snakemake workflow with STAR and RSEM. It can be found [here](https://github.com/tjbencomo/nmsc-star).

## Data
### Processed Data
The processed data for this project is available at [Zenodo](https://zenodo.org/records/10272679).

### Raw Data
Raw data can easily be downloaded from the apropriate sites with the data download scripts:

* `download_sra.sh` - download all FASTQ files available on SRA. Pulls SRR IDs from `srr_ids.txt`
* `download_ena.sh` - download all FASTQ files available on BioStudies. Pulls ERR IDs from `ena_fastq_urls.txt`

`srr_ids.txt` and `ena_fastq_urls.txt` can be found on Zenodo. Once the FASTQ files are downloaded, the Snakemake pipeline mentioned
above will generate BAM files and DESeq2 analysis

## Analysis Scripts


