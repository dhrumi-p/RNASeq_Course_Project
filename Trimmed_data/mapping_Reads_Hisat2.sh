#!/bin/bash

#SBATCH --job-name=mapping_hisat2_trimmed_files
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=30G
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/mapping_hisat2_trimmed_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/mapping_hisat2_trimmed_%j.e
#SBATCH --mail-user=dhrumi.patel@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


# Defining paths as variables 
# path to container
container="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
# samtools 1.20 (version)
reference_dir="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome"
# below is location to where my index file for the mus musculus is stored 
hisat2_index="$reference_dir/Mus_musculus_DNA_index"
# reads_dir is the path where my fastq files are stored and output_dir is the path to where i want to store the output  results 
reads_dir="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed"   
output_dir="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/mapped_files_trimmed"


# Loop over all forward reads (ending with _1.fastq.gz)
for R1 in "$reads_dir"/*_1.trimmed.fastq.gz
do
    # Extract sample name
    sample=$(basename "$R1" _1.trimmed.fastq.gz)
    
    # Define reverse read and output SAM file and for output .sam extension i have taken bcz it creates the SAM files 
    R2="$reads_dir/${sample}_2.trimmed.fastq.gz"
    OUT="$output_dir/${sample}.trimmed.sam"

    # Aligining the pair-end RNA-seq reads to indexed reference genome using, Hisat2
    # i have specified the strandness as RF (reverse-forward) since i have pair-end data 
    apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples "$container" \
        hisat2 -p $SLURM_CPUS_PER_TASK \
        --rna-strandness RF \
        -x "$hisat2_index" \
        -1 "$R1" -2 "$R2" \
        -S "$OUT"
done


container_sam="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# conversion of SAM file to BAM files 
for SAM_FILE in "$output_dir"/*.trimmed.sam
do
    sample=$(basename "$SAM_FILE" .trimmed.sam) # getting the sample name 
    BAM_FILE="$output_dir/${sample}.trimmed.bam" # giving the path where BAM files needs to be stored 
    
    # Convert SAM to BAM using samtools 
    apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples "$container_sam" \
        samtools view -bS "$SAM_FILE" > "$BAM_FILE"

done











