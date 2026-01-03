#!/bin/bash
#SBATCH --job-name=adapter_trimming
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/trimming_Counts_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/trimming_Counts_%j.e

# defining paths   
CONTAINER="/containers/apptainer/fastp_0.24.1.sif"
#  fastp v0.24.1
RAW_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/reads_Blood"  # folder containing *.fastq.gz (sample counts file)
OUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed" # folder for trimmed files


# processing each sample   
for R1 in "$RAW_DIR"/*_1.fastq.gz; do 
    
    # Derive R2 filename
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    # Extract sample name (removes path + R1)
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    # Output file names 
    OUT_R1="$OUT_DIR/${SAMPLE}_1.trimmed.fastq.gz"
    OUT_R2="$OUT_DIR/${SAMPLE}_2.trimmed.fastq.gz"
    HTML="$OUT_DIR/${SAMPLE}_fastp.html"
    JSON="$OUT_DIR/${SAMPLE}_fastp.json"

# To do quality control and adapter trimming on pair-end RNAseq reads
    apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples "$CONTAINER" fastp \
        -i "$R1" -I "$R2" \
        -o "$OUT_R1" -O "$OUT_R2" \
        -w $SLURM_CPUS_PER_TASK \
        --detect_adapter_for_pe \
        --thread $SLURM_CPUS_PER_TASK \
        --html "$HTML" \
        --json "$JSON"

done

# last two lines is to generate the HTML and JSON QC report 