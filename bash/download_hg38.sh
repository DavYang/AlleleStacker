#!/bin/bash
#SBATCH --job-name=hg38_download
#SBATCH --time=2-00:00:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=download.%A_%a.out


mkidir -p /reference_genomes/hg38
cd /reference_genomes/hg38

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"

wget "${BASE_URL}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
wget "${BASE_URL}/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" 

echo "done"
