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

# Function to parse YAML
parse_yaml() {
    local yaml_file=$1
    local prefix=$2
    local s
    local w
    local fs

    s='[[:space:]]*'
    w='[a-zA-Z0-9_.-]*'
    fs=$(echo @|tr @ '\034')

    (
        sed -e '/- [^\"]'"[^\']"'.*: /s|\([ ]*\)- \([[:space:]]*\)|\1-\'$'\n''  \1\2|g' |
        sed -ne '/^--/s|--||g; s|\"|\\\"|g; s/[[:space:]]*$//g;' \
            -e 's/\$/\\\$/g' \
            -e "/#.*"/!s| #.*$||g; s|^\($s\)\($w\)$s:$s\"\(.*\)\"$|\1$fs\2$fs\3|p" \
            -e "s|^\($s\)\($w\)${s}[:-]$s\(.*\)$|\1$fs\2$fs\3|p" |
        awk -F$fs '{
            indent = length($1)/2;
            if (length($2) == 0) { conj[indent]="+";} else {conj[indent]="";}
            vname[indent] = $2;
            for (i in vname) {if (i > indent) {delete vname[i]}}
            if (length($3) > 0) {
                vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
                printf("%s%s%s%s=(\"%s\")\n", "'"$prefix"'",vn, $2, conj[indent-1], $3);
            }
        }'
    ) < "$yaml_file"
}

# Parse config and create variables
eval $(parse_yaml "$CONFIG_FILE")

# Validate required config values
required_vars=("reference_fasta" "python" "vcf_base_dir" "bed_base_dir" "conda_analysis_env")
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
