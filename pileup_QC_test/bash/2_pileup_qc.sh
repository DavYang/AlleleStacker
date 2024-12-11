#!/bin/bash
#SBATCH --partition=unlimited
#SBATCH --time=infinite
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=20
#SBATCH --output=logs.qc_%j.out
#SBATCH --error=logs.qc_%j.err

# Command line arguments
SAMPLE=$1
VCF_DIR=$2
REF_FASTA=$3
BED_DIR=$4
OUTPUT_DIR=$5



mkdir -p ${OUTPUT_DIR}
# Load conda
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate bcftools


# conda deactivate
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Run QC script
python python/12-10_test.py \
    --vcf ${VCF_DIR}/${SAMPLE}.region.vcf.gz \
    --ref ${REF_FASTA} \
    --prefix ${SAMPLE} \
    --bed-dir ${BED_DIR} \
    --outdir ${OUTPUT_DIR} \
    --threads 20