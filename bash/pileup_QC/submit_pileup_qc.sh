#!/bin/bash
# submit_pileup_qc.sh

SAMPLES_FILE=$1
CONFIG_FILE=$2
BASE_OUTPUT_DIR=$3
PYTHON_DIR=$4

# Check args
if [ $# -ne 4 ]; then
    echo "Usage: $0 SAMPLES_FILE CONFIG_FILE OUTPUT_DIR"
    echo "Example: $0 samples.txt config.yaml /path/to/output"
    exit 1
fi

# Verify inputs
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE not found"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file $CONFIG_FILE not found"
    exit 1
fi

# Parse config file
while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*#.*$ ]] && continue
    
    # Remove leading/trailing whitespace and comments
    line=$(echo "$line" | sed 's/#.*$//' | xargs)
    [[ -z "$line" ]] && continue
    
    # Check if line contains a colon (key-value pair)
    if [[ "$line" == *":"* ]]; then
        # Extract key and value
        key=$(echo "$line" | cut -d: -f1 | xargs)
        value=$(echo "$line" | cut -d: -f2- | sed 's/^[[:space:]]*"//' | sed 's/"[[:space:]]*$//' | xargs)
        
        # Convert key format (remove spaces and dashes)
        key=$(echo "$key" | tr -d ' ' | tr '-' '_')
        
        # Special handling for nested conda settings
        if [[ "$key" == "conda" ]]; then
            in_conda_section=true
            continue
        fi
        
        if [[ "$in_conda_section" == true ]]; then
            if [[ "$key" == "analysis_env" ]]; then
                declare "conda_analysis_env=$value"
                in_conda_section=false
            fi
        else
            # Regular variable declaration
            declare "$key=$value"
        fi
    fi
done < "$CONFIG_FILE"

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
        --error=:"logs/${SAMPLE}_%j.err" \
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