#!/bin/bash
#SBATCH --job-name=multiqc_qc
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/multiqc_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/multiqc_%j.e


# Defining paths

CONTAINER="/containers/apptainer/multiqc-1.19.sif"
#multiqc v - 1.19

# Directory containing FastQC / fastp outputs
RAW_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed"

# Output directory for MultiQC report
OUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed/Multiqc"


# Multi-qc
apptainer exec \
    --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples \
    "$CONTAINER" \
    multiqc "$RAW_DIR" \
    --outdir "$OUT_DIR" \
    --filename multiqc_report.html \
    
