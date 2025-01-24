#!/bin/bash                                                                   
#SBATCH --job-name=plot_unmethylated_dist_sample                              
#SBATCH --output=logs/plot_unmethylated_dist_sample_%A_%a.out                 
#SBATCH --error=logs/plot_unmethylated_dist_sample_%A_%a.err                  
#SBATCH --array=1-21 # Adjust as needed                                       
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
                                                                              
# Get sorted list of sample directories, excluding specific subdirectories    
cd "$1"

SAMPLES=($(ls -d "$1"/SPM* | grep -v 
"log_files\|regions_by_label\|segmentation_scripts\|unmethylated_distribution_
plots" | xargs -n1 basename | sort))
CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))                                    
SAMPLE=${SAMPLES[$CURRENT_INDEX]}                                             
                                                                              
# Path to the Python script                                                   
PYTHON="./python/segmentation"                                                
                                                                              
# Run the Python plotting script                                              
python "$PYTHON/plot_unmethylated_distribution.py" "$1" "$2" "$SAMPLE"  
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
    --input-dir "$1/regions" \
    --output-dir "$2" \
    --sample-name "$SAMPLE"
