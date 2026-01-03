#!/bin/bash

#SBATCH --job-name=sort_Trim_BAM_Trimmed
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=25GB
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/Sort_Trim_BAM_Trimmed_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/Sort_Trim_BAM_Trimmed_%j.e
#SBATCH --mail-user=dhrumi.patel@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

# Defining paths 
INPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/mapped_files_trimmed"
CONTAINER_SAM="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
# samtools 1.20 (version)
OUTPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/mapped_files_trimmed"


# using for loop to iterate throught all the BAM files present in the input folder 
for bam in "$INPUT_DIR"/*.trimmed.bam; do
    base=$(basename "$bam" .trimmed.bam) # getting the name of file , which later to be used for the naming the output files
    sorted_bam="$OUTPUT_DIR/${base}_sorted.trimmed.bam" # defining the output file name

    apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples \
    "$CONTAINER_SAM" samtools sort -o "$sorted_bam" "$bam" # sorting BAM files by genomic coordinates using samtools

    apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples \
    "$CONTAINER_SAM" samtools index "$sorted_bam" # indexing the BAM file using samtools 

done



