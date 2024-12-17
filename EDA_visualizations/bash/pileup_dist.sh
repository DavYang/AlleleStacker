#!/bin/bash
#SBATCH --job-name=pileup_dist
#SBATCH --output=logs/pileup_dist_%A_%a.out
#SBATCH --error=logs/pileup_dist_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Set the base directory
BASE_DIR="$1"

# Set the output base directory
OUTPUT_BASE_DIR="$2"

# Add debugging information
echo "Current working directory: $(pwd)"
echo "Base directory: ${BASE_DIR}"
echo "Output base directory: ${OUTPUT_BASE_DIR}"

# Get the list of sample directories
SAMPLE_DIRS=($(ls -d ${BASE_DIR}/SPM* 2>/dev/null))


# Get the current sample directory based on the array task ID
SAMPLE_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID-1]}"

# Add debugging information
echo "Current sample directory: ${SAMPLE_DIR}"

# Set the output directory
OUTPUT_DIR="${OUTPUT_BASE_DIR}/$(basename ${SAMPLE_DIR})"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Add debugging information
echo "Output directory: ${OUTPUT_DIR}"

# Run the Python script
echo "Running Python script..."
python ./python/pileup_dist.py ${SAMPLE_DIR} ${OUTPUT_DIR}
