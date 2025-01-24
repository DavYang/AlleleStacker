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
