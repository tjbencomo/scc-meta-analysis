#!/bin/bash
#SBATCH --job-name=star-fusion
#SBATCH --output=/scratch/users/tbencomo/scc-cell-lines/log
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu


set -e

cd $SCRATCH/scc-cell-lines

star_env=$CONTAINERS/star-fusion.v1.12.0.simg
ref_dir=$GROUP_HOME/refs/rnaseq-refs/ctat/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/

echo "Running A431-SRR17247482"
singularity exec $star_env STAR-Fusion \
    --genome_lib_dir $ref_dir \
    --left_fq fastq/SRR17247482_1.fastq \
    --right_fq fastq/SRR17247482_2.fastq \
    --output_dir A431-SRR17247482 \
    --STAR_SortedByCoordinate \
    --STAR_outSAMattrRGline "ID:id SM:A431-SRR17247482" \
    --examine_coding_effect \
    --FusionInspector validate \
    --CPU 8

echo "Running SCCIC1-SRR8286261"
singularity exec $star_env STAR-Fusion \
    --genome_lib_dir $ref_dir \
    --left_fq fastq/SRR8286261_1.fastq \
    --right_fq fastq/SRR8286261_2.fastq \
    --output_dir SCCIC1-SRR8286261 \
    --STAR_SortedByCoordinate \
    --STAR_outSAMattrRGline "ID:id SM:SCCIC1-SRR8286261" \
    --examine_coding_effect \
    --FusionInspector validate \
    --CPU 8

echo "Running SCCIC1-SRR8286262"
singularity exec $star_env STAR-Fusion \
    --genome_lib_dir $ref_dir \
    --left_fq fastq/SRR8286262_1.fastq \
    --right_fq fastq/SRR8286262_2.fastq \
    --output_dir SCCIC1-SRR8286262 \
    --STAR_SortedByCoordinate \
    --STAR_outSAMattrRGline "ID:id SM:SCCIC1-SRR8286262" \
    --examine_coding_effect \
    --FusionInspector validate \
    --CPU 8


echo "Running A431-SRR12453889"
# singularity exec $star_env STAR-Fusion \
#     --genome_lib_dir $ref_dir \
#     --left_fq fastq/SRR12453889_1.fastq \
#     --right_fq fastq/SRR12453889_2.fastq \
#     --output_dir A431-SRR12453889 \
#     --STAR_SortedByCoordinate \
#     --STAR_outSAMattrRGline "ID:id SM:A431-SRR12453889" \
#     --examine_coding_effect \
#     --FusionInspector validate \
#     --CPU 8

echo "Running A431-SRR12453888"
# singularity exec $star_env STAR-Fusion \
#     --genome_lib_dir $ref_dir \
#     --left_fq fastq/SRR12453888_1.fastq \
#     --right_fq fastq/SRR12453888_2.fastq \
#     --output_dir A431-SRR12453888 \
#     --STAR_SortedByCoordinate \
#     --STAR_outSAMattrRGline "ID:id SM:A431-SRR12453888" \
#     --examine_coding_effect \
#     --FusionInspector validate \
#     --CPU 8

