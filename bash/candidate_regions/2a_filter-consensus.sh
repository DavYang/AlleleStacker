#!/bin/bash
#SBATCH --job-name=filter-consensus
#SBATCH --time=4:00:00
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5  
#SBATCH --mem=64gb         
#SBATCH --output=filter-consensus_%A.out
#SBATCH --error=filter-consensus_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.yang2@einsteinmed.edu

# Enable CPU binding for better thread performance
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Add explicit CPU affinity for Python processes
export PYTHONPATH="${PYTHONPATH}:$PYTHON"
export NUMEXPR_MAX_THREADS=5
export MKL_NUM_THREADS=5
export OPENBLAS_NUM_THREADS=5

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Set script parameters

REGIONS_DIR="$1"
INPUT_DIR="$2"
OUTPUT_DIR="$3"
# Get haplotype from command line argument
HAPLOTYPE="$4"

SAMPLE_LIST="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/all_samples.txt"
PYTHON="./python"



mkdir -p $OUTPUT_DIR



if [ -z "$HAPLOTYPE" ]; then
    echo "Error: Haplotype not specified"
    echo "Usage: sbatch script.sh H1|H2"
    exit 1
fi

echo "Starting processing for $HAPLOTYPE with $(nproc) available CPUs"

# Run the Python script
python $PYTHON/filter_consensus.py \
    --consensus_regions "$INPUT_DIR/consensus_regions_${HAPLOTYPE}_consensus.bed" \
    --regions_dir "$REGIONS_DIR/${HAPLOTYPE}_M" \
    --output_file "$OUTPUT_DIR/filtered_consensus_${HAPLOTYPE}.bed" \
    --sample-list-file "$SAMPLE_LIST" \
    --haplotype "$HAPLOTYPE" \
    --min-samples 1 #adjust as needed

echo "Completed processing for $HAPLOTYPE"

