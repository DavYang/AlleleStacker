#!/bin/bash
#SBATCH --job-name=submit-jobs
#SBATCH --time=00:10:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --output=submit_jobs.%A_%a.out

# Define the paths
main_dir="/gs/gsfs0/users/greally-lab/David/6_base-seq_SC1/WGS-analysis/outputs_compiled/20240501-results/cpg_pileup_beds"
output_dir="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/outputs/segmentation_regions"
log_out="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/outputs/segmentation_regions/log_files"
scripts_dir="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/outputs/segmentation_regions/segmentation_scripts"

mkdir -p "$output_dir"
mkdir -p "$log_out"
mkdir -p "$scripts_dir"

methbat="/gs/gsfs0/shared-lab/greally-lab/David/software/methbat-v0.13.2-x86_64-unknown-linux-gnu/methbat"

# Iterate through each sample directory
for sample_dir in "$main_dir"/*; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi

    # Get sample name from the directory name
    sample_name=$(basename "$sample_dir")

    # Create a unique job name
    job_name="segment-${sample_name}"

    # Create a Slurm script for each sample
    slurm_script="$scripts_dir/${job_name}.sh"

    cat > "$slurm_script" <<EOL
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --time=2-00:00:00 
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=${log_out}/${job_name}.segment.%A_%a.out

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate DYenv

# Change to the sample directory
cd "${sample_dir}"

# Run segmentation
bed_files=(*.bed)

# Check if there are any bed files
if [ \${#bed_files[@]} -eq 0 ]; then
    echo "No bed files found in the sample directory."
    exit 1
fi

echo "Running segmentation for \${bed_files[@]}"
"${methbat}" segment --input-prefix "${sample_name}.GRCh38" --output-prefix "${output_dir}/${sample_name}" --enable-haplotype-segmentation --condense-bed-labels --min-cpgs 5 --max-gap 500

if [ \$? -ne 0 ]; then
    echo "Error: command failed for \${bed_files[@]} with exit code \$?"
else
    echo "Segmentation completed successfully for \${bed_files[@]}"
fi
EOL

    # Submit the Slurm job for each sample
    sbatch "$slurm_script"
done
