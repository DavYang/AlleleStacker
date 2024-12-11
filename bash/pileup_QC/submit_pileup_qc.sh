#!/bin/bash
# submit_pileup_qc.sh

SAMPLES_FILE=$1
CONFIG_FILE=$2
BASE_OUTPUT_DIR=$3

# Check args
if [ $# -ne 3 ]; then
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

# Load config - Fixed version
eval "$(python3 -c 'import yaml;
import json;
with open("'"$CONFIG_FILE"'") as f:
    cfg = yaml.safe_load(f);
for k,v in cfg.items():
    if isinstance(v, dict):
        for sk,sv in v.items():
            print(f"{k}_{sk}={json.dumps(str(sv))}");
    else:
        print(f"{k}={json.dumps(str(v))}");
')" || {
    echo "Error loading config file"
    exit 1
}

# Create base directories
mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${BASE_OUTPUT_DIR}/logs"

# Submit jobs
while read SAMPLE; do
    [[ -z "${SAMPLE}" || "${SAMPLE}" =~ ^#.*$ ]] && continue
    
    # Define sample directory before submission
    SAMPLE_DIR="${BASE_OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_DIR}"


    sbatch \
        --job-name="${SAMPLE}" \
        --partition=normal \
        --time=2-00:00:00 \
        --mem=16G \
        --cpus-per-task=20 \
        --output="${BASE_OUTPUT_DIR}/logs/${SAMPLE}_%j.out" \
        --error="${BASE_OUTPUT_DIR}/logs/${SAMPLE}_%j.err" \
        --wrap="
            set -e
            source /public/apps/conda3/etc/profile.d/conda.sh

            # Create sample-specific directories
            SAMPLE_DIR=${BASE_OUTPUT_DIR}/${SAMPLE}

            # Run QC
            conda activate ${conda_analysis_env}
            python python/pileup_QC.py \
                --vcf ${vcf_base_dir}/${SAMPLE}/${SAMPLE}.*.vcf.gz \
                --ref ${reference_fasta} \
                --prefix ${SAMPLE} \
                --bed-dir ${bed_base_dir}/${SAMPLE} \
                --outdir ${SAMPLE_DIR} \
                --threads 20
            "
    
    echo "Submitted job for ${SAMPLE} (output in ${SAMPLE_DIR})"
done < "${SAMPLES_FILE}"