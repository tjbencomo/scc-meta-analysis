#!/bin/bash
#SBATCH --job-name=download_ena
#SBATCH --output=
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

set -euo pipefail

ml system parallel
# fqdir=$SCRATCH/nmsc-data
fqdir=$GROUP_SCRATCH/nmsc-data
mkdir -p $fqdir
cd $fqdir
ena_file=$HOME/nmsc-rna-seq/data/ena_fastq_urls.txt
cat $ena_file | parallel -j 24 'wget {}'
