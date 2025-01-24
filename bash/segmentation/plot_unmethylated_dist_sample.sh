#!/bin/bash
#SBATCH --job-name=unmeth_dist
#SBATCH --output=logs/unmeth_dist_%A_%a.out
#SBATCH --error=logs/unmeth_dist_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Create output and log directories
mkdir -p "$2"
mkdir -p logs

# Get list of unique sample names from the bed files in H1_M directory
# Using proper path construction and error handling
input_dir="$1/regions/H1_M"
if [ ! -d "$input_dir" ]; then
    echo "Error: Directory not found: $input_dir"
    exit 1
fi

cd "$input_dir" || exit 1

# Get sample names from bed files, removing the _H1_M.bed suffix
SAMPLES=($(ls -1 SPM*.bed 2>/dev/null | sed 's/_H1_M\.bed$//' | sort -u))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: No SPM*.bed files found in $input_dir"
    exit 1
fi

CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [ $CURRENT_INDEX -ge ${#SAMPLES[@]} ]; then
    echo "Error: Array index $CURRENT_INDEX is out of bounds. Only ${#SAMPLES[@]} samples found."
    exit 1
fi

SAMPLE=${SAMPLES[$CURRENT_INDEX]}

echo "Processing sample: $SAMPLE"
echo "Input directory: $1"
echo "Output directory: $2"

# Run the Python script with the correct paths
python ./python/segmentation/plot_unmethylated_distribution.py \
    --input-dir "$1/regions" \
    --output-dir "$2" \
    --sample-name "$SAMPLE"
#!/bin/bash
#SBATCH --job-name=unmeth_dist
#SBATCH --output=logs/unmeth_dist_%A_%a.out
#SBATCH --error=logs/unmeth_dist_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Create output and log directories
mkdir -p "$2"
mkdir -p logs

# Get list of unique sample names from the bed files in H1_M directory
# This assumes all samples have files in H1_M
SAMPLES=($(ls "$1"/regions/H1_M/SPM*.bed | xargs -n1 basename | sed 's/_H1_M.bed//g' | sort -u))
CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$CURRENT_INDEX]}

# Run the Python script with the correct paths
python ./python/segmentation/plot_unmethylated_distribution.py \
    --input-dir "$1" \
    --output-dir "$2" \
    --sample-name "$SAMPLE"
