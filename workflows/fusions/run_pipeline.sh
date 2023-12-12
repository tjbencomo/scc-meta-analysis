#!/bin/bash
#SBATCH --job-name=nmsc-fusions
#SBATCH --output=/home/users/tbencomo/nmsc-rna-seq/workflows/fusions/log
#SBATCH --nodes=1
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu

# OUT OF USE: SBATCH --qos=long

set -e
cd /home/users/tbencomo/nmsc-rna-seq/workflows/fusions
snakemake --cluster-config cluster.json -j 499 \
    --rerun-incomplete \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
