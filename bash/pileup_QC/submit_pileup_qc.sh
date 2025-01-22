#!/bin/bash
# submit_pileup_qc.sh

SAMPLES_FILE=$1
CONFIG_FILE=$2
BASE_OUTPUT_DIR=$3
PYTHON_DIR=$4

# Check args
if [ $# -ne 4 ]; then
    echo "Usage: $0 SAMPLES_FILE CONFIG_FILE OUTPUT_DIR PYTHON_DIR"
    echo "Example: $0 samples.txt config.yaml output_dir python/"
    exit 1
fi

# Convert to absolute paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SAMPLES_FILE="$(realpath "$1")"
CONFIG_FILE="$(realpath "$2")"
BASE_OUTPUT_DIR="$(realpath "$3")"
PYTHON_DIR="$(realpath "$4")"

# Verify inputs
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE not found"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found"
    exit 1
fi

# Function to parse YAML
parse_yaml() {
    local yaml_file=$1
    local prefix=$2
    while IFS=':' read -r key value; do
        # Skip comments and empty lines
        [[ $key =~ ^[[:space:]]*# ]] && continue
        [[ -z "${key}" ]] && continue
        
        # Remove leading/trailing whitespace and comments
        key=$(echo "$key" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        value=$(echo "$value" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sed 's/[[:space:]]*#.*$//')
        
        # Skip empty values
        [[ -z "${value}" ]] && continue
        
        # Handle nested values under conda:
        if [[ "$key" == "conda" ]]; then
            in_conda=true
            continue
        fi
        if [[ "$in_conda" == true ]]; then
            if [[ "$key" =~ ^[[:space:]] ]]; then
                # Remove leading spaces and combine with conda prefix
                key="conda_$(echo "$key" | sed 's/^[[:space:]]*//')"
            else
                in_conda=false
            fi
        fi
        
        # Export the variable
        echo "${prefix}${key}=\"${value}\""
    done < "$yaml_file"
}

# Parse config and create variables
eval $(parse_yaml "$CONFIG_FILE")

# Validate required config values
required_vars=("reference_fasta" "vcf_base_dir" "bed_base_dir" "conda_analysis_env")
for var in "${required_vars[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Required config variable $var is not set"
        exit 1
    fi
done

# Create base directories
mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p logs

# Submit jobs
while read SAMPLE; do
    [[ -z "${SAMPLE}" || "${SAMPLE}" =~ ^#.*$ ]] && continue
    
    # Define sample directory before submission
    SAMPLE_DIR="${BASE_OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_DIR}"

    sbatch \
        --job-name="${SAMPLE}" \
        --partition=quick \
        --time=4:00:00 \
        --mem=32G \
        --cpus-per-task=40 \
        --output="logs/${SAMPLE}_%j.out" \
        --error="logs/${SAMPLE}_%j.err" \
        --wrap="
            set -e

            # Create sample-specific directories
            SAMPLE_DIR=${BASE_OUTPUT_DIR}/${SAMPLE}

            # Run QC
            source /public/apps/conda3/etc/profile.d/conda.sh
            conda activate ${conda_analysis_env}
            python ${PYTHON_DIR}/pileup_QC.py \
                --vcf ${vcf_base_dir}/${SAMPLE}/${SAMPLE}.*.vcf.gz \
                --ref ${reference_fasta} \
                --prefix ${SAMPLE} \
                --bed-dir ${bed_base_dir}/${SAMPLE} \
                --outdir ${SAMPLE_DIR} \
                --threads 40
            "
    
    echo "Submitted job for ${SAMPLE} (output in ${SAMPLE_DIR})"
done < "${SAMPLES_FILE}"
