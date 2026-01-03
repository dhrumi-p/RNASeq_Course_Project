#!/bin/bash

#SBATCH --job-name=cleaning_counts_data
#SBATCH --cpus-per-task=6
#SBATCH --time=01:00:00
#SBATCH --mem=6GB
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/scripts_blood/output/Cleaning_Counts_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/scripts_blood/error/Cleaning_Counts_%j.e


# Input file
INPUT=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/gene_counts_trimmed.txt
STEP1=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/trimmed_gene_counts_step1.txt
STEP2=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/trimmed_gene_counts_step2.txt
STEP3=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/trimmed_gene_counts_step3.txt
OUTPUT=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/gene_counts_trimmed.csv

# here, i am preparing the counts_data for the downstream processing 
# Step 1: Remove comment lines that is the first line from the feature counts file 
grep -v "^#" "$INPUT" > "$STEP1"

# Step 2: Remove annotation columns (keep Geneid(1st column) + counts(columns from 8 onwards)) (remove col2 to 7)
cut -f1,7- "$STEP1" > "$STEP2"

# Step 3: Clean sample names (remove paths and .bam) (cleaning the column names from 8 th column, just to have sampel names, thats all)
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=2;i<=NF;i++){sub(".*/","",$i); sub("_sorted.bam","",$i)}}{print}' "$STEP2" > "$STEP3"

# Step 4: Convert tab separated file(tsv) to CSV (comma separated file)
sed 's/\t/,/g' "$STEP3" > "$OUTPUT"
