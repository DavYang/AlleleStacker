#!/bin/bash
#SBATCH --job-name=pile
#SBATCH --time=4:00:00
#SBATCH --partition=quick
#SBATCH --mail-type=END
#SBATCH --mail-user=david.yang2@einsteinmed.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=log.pile_%j.out
#SBATCH --error=log.pile_%j.err

# create conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate bcftools

# Command line arguments
VCF_FILE=$1
BAM_FILE=$2
REF_FASTA=$3
OUTPUT_PREFIX=$4
OUTPUT_DIR=$5
CPG_TOOLS_DIR=$6  

# Create output directories
mkdir -p ${OUTPUT_DIR}

# Step 1: Extract test region
REGION="chr1:167535318-167685531"


echo "Extracting region ${REGION}..."

# Extract and index region from VCF
bcftools view -r ${REGION} ${VCF_FILE} -Oz -o ${OUTPUT_DIR}/${OUTPUT_PREFIX}.region.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${OUTPUT_PREFIX}.region.vcf.gz


conda deactivate
conda activate biobase
# Extract and index region from BAM
samtools view -b ${BAM_FILE} ${REGION} > ${OUTPUT_DIR}/${OUTPUT_PREFIX}.region.bam
samtools index ${OUTPUT_DIR}/${OUTPUT_PREFIX}.region.bam

# Step 2: Run pb-CpG-tools v2.3.2 (default settings listed out)
echo "Running pb-CpG-tools v2.3.2..."

${CPG_TOOLS_DIR}/bin/aligned_bam_to_cpg_scores \
    --bam ${OUTPUT_DIR}/${OUTPUT_PREFIX}.region.bam \
    --pileup-mode model \
    --modsites-mode denovo \
    --model ${CPG_TOOLS_DIR}/models/pileup_calling_model.v1.tflite \
    --ref ${REF_FASTA} \
    --output-prefix ${OUTPUT_DIR}/${OUTPUT_PREFIX} \
    --threads 8 \
    --min-coverage 4 \
    --min-mapq 1 \
    --hap-tag HP 

echo "CpG analysis complete."
