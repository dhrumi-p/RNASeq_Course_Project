#!/bin/bash

#SBATCH --job-name=Quality_control_trimmed_data
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/QC_trimmed_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/QC_trimmed_%j.e
#SBATCH --mail-user=dhrumi.patel@students.unibe.ch
#SBATCH --mail-type=begin,end,fail 

 # setting paths, to acess container and also the reads file and where to get the output 
 CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"
 # fastqc v0.12.1
 INPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed"
 OUTPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/QC_Trimmed"

 # to run containers apptainers are required 
 # By default apptainer reads my directories to the container, it does for home directory, but that doesnto have my sample
 # and the other folders and so bind is used here 
 # so binds helps to connect the folder from my computer(here cluster) to the container
 # so, first we have added the input directory in code and then with -o we have given the path for output directory

# Here, i am trying to do the QC on the already trimmed file to visualize the difference 

 apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples "$CONTAINER" fastqc \
    "$INPUT_DIR"/*.trimmed.fastq.gz \
     -o "$OUTPUT_DIR" \
     -t $SLURM_CPUS_PER_TASK 