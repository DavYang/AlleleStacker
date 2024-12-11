#!/bin/bash
#SBATCH --job-name=segment
#SBATCH --time=4:00:00
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --output=log.segment.%A_%a.out
#SBATCH --error=log.segment.%A_%a.err

# Define the paths
input_dir="$1"
sample_name="$2"
output_dir="$3"

mkdir -p "${output_dir}"

methbat="/gs/gsfs0/shared-lab/greally-lab/David/software/methbat-v0.13.2-x86_64-unknown-linux-gnu/methbat"

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# "${methbat}" segment --input-prefix "${input_dir}/${sample_name}" --output-prefix "${output_dir}/${sample_name}" --enable-haplotype-segmentation --condense-bed-labels --min-cpgs 2 --max-gap 200

"${methbat}" segment --input-prefix "${input_dir}/${sample_name}" --output-prefix "${output_dir}/${sample_name}" --enable-haplotype-segmentation --condense-bed-labels --min-cpgs 2 --max-gap 200