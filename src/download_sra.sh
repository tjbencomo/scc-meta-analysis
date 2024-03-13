#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=
#SBATCH --nodes=1
#SBATCH --time=02-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=14000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=

set -euo pipefail

ml system parallel
fqdir=$SCRATCH/nmsc-data
mkdir -p $fqdir
cd $fqdir
srr_file=$HOME/scc-meta-analysis/data/srr_ids.txt
cat $srr_file | parallel -j 12 'fasterq-dump {}'
