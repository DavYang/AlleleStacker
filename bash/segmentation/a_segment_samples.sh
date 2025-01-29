#!/bin/bash
#SBATCH --job-name=submit-jobs
#SBATCH --time=00:10:00
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --output=submit_jobs.%A_%a.out
#SBATCH --error=submit_jobs.%A_%a.err

# Enable debug output
set -x

# Accept parameters
MAIN_DIR="$1"
OUTPUT_DIR="$2"
LOG_OUT="${OUTPUT_DIR}/log_files"
SCRIPTS_DIR="${OUTPUT_DIR}/segmentation_scripts"

# Validate required parameters
if [ -z "$MAIN_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Both input and output directories must be specified"
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_OUT"
mkdir -p "$SCRIPTS_DIR"

# Print debug info
echo "Main directory: $MAIN_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOG_OUT"
echo "Scripts directory: $SCRIPTS_DIR"

# Path to methbat
methbat="/gs/gsfs0/shared-lab/greally-lab/David/software/methbat-v0.13.2-x86_64-unknown-linux-gnu/methbat"

# Iterate through each sample directory
cd "$MAIN_DIR" || exit 1
sample_dirs=$(find . -maxdepth 1 -type d -name 'SPM*')
echo "Found sample directories: $sample_dirs"

if [ -z "$sample_dirs" ]; then
    echo "No SPM directories found"
    exit 0
fi

for sample_dir in $sample_dirs; do
    sample_dir="${MAIN_DIR}/${sample_dir#./}"
    echo "Processing directory: $sample_dir"

    # Get sample name from the directory name
    sample_name=$(basename "$sample_dir")
    echo "Sample name: $sample_name"
    
    # Check if passing directory exists
    passing_dir="${sample_dir}/filtered_beds/passing"
    echo "Checking passing directory: $passing_dir"
    
    if [ ! -d "$passing_dir" ]; then
        echo "Warning: No passing directory found for ${sample_name}"
        continue
    fi

    # Create sample-specific output directory
    sample_output_dir="${OUTPUT_DIR}/${sample_name}"
    mkdir -p "$sample_output_dir"
    echo "Created output directory: $sample_output_dir"

    # Create a unique job name
    job_name="${sample_name}-segment"
    echo "Job name: $job_name"

    # Create a Slurm script for each sample
    slurm_script="$SCRIPTS_DIR/${job_name}.sh"
    echo "Creating script: $slurm_script"

    cat > "$slurm_script" <<'EOL'
#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --time=4:00:00 
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --output=LOGOUT/JOBNAME.segment.%A_%a.out

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Create temporary directory for renamed files
temp_dir="OUTPUTDIR/SAMPLENAME"
mkdir -p "$temp_dir"

# Process haplotype 1
if [ -f "PASSINGDIR/SAMPLENAME.filtered.hap1.bed" ]; then
    # Create renamed copy for haplotype 1
    cp "PASSINGDIR/SAMPLENAME.filtered.hap1.bed" "$temp_dir/SAMPLENAME.haplotype1.combined.bed"
    
    echo "Running segmentation for SAMPLENAME haplotype 1"
    cd "$temp_dir"
    METHBAT segment \
        --input-prefix "SAMPLENAME.haplotype1" \
        --output-prefix "OUTPUTDIR/SAMPLENAME.hap1" \
        --condense-bed-labels \
        --enable-nodata-segments \
        --min-cpgs 5 \
        --max-gap 500
    
    if [ $? -ne 0 ]; then
        echo "Error: Segmentation failed for SAMPLENAME haplotype 1"
    fi
fi

# Process haplotype 2
if [ -f "PASSINGDIR/SAMPLENAME.filtered.hap2.bed" ]; then
    # Create renamed copy for haplotype 2
    cp "PASSINGDIR/SAMPLENAME.filtered.hap2.bed" "$temp_dir/SAMPLENAME.haplotype2.combined.bed"
    
    echo "Running segmentation for SAMPLENAME haplotype 2"
    cd "$temp_dir"
    METHBAT segment \
        --input-prefix "SAMPLENAME.haplotype2" \
        --output-prefix "OUTPUTDIR/SAMPLENAME.hap2" \
        --enable-nodata-segments \
        --condense-bed-labels \
        --min-cpgs 5 \
        --max-gap 500
    
    if [ $? -ne 0 ]; then
        echo "Error: Segmentation failed for SAMPLENAME haplotype 2"
    fi
fi

# Cleanup temporary files
rm -rf "$temp_dir"

echo "Processing completed for SAMPLENAME"
EOL

    # Replace placeholders with actual values
    sed -i "s|JOBNAME|${job_name}|g" "$slurm_script"
    sed -i "s|LOGOUT|${LOG_OUT}|g" "$slurm_script"
    sed -i "s|SAMPLENAME|${sample_name}|g" "$slurm_script"
    sed -i "s|PASSINGDIR|${passing_dir}|g" "$slurm_script"
    sed -i "s|OUTPUTDIR|${sample_output_dir}|g" "$slurm_script"
    sed -i "s|METHBAT|${methbat}|g" "$slurm_script"

    # Make script executable
    chmod +x "$slurm_script"
    echo "Created script: $slurm_script"

    # Submit the Slurm job for each sample
    sbatch "$slurm_script"
done

echo "Script generation complete"