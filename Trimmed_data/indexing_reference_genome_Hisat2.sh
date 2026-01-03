#!/bin/bash

#SBATCH --job-name=indexing_hisat2
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem=8000M
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/scripts_blood/output/indexing_hisat2_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/scripts_blood/error/indexing_hisat2_%j.e
#SBATCH --mail-user=dhrumi.patel@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


# setting the paths to the container and wher the files can be found respectigely how the index files should be named.
container="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
# hisat2 2.2.1 (version of hisat2)
reference_dir="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome"
genome_fasta="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa"
gtf_file="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome/Mus_musculus.GRCm39.115.gtf"
index_prefix="$reference_dir/Mus_musculus_DNA_index"

# creating splice site and exon files from the reference genome using annotation(.gtf) file 
# after >(redirection operator) we have given the file name and a path ahead of it where exactly the new file
# namely Mus_musculus_splicesites.txt and Mus_musculus_exons.txt be created 

apptainer exec --bind $reference_dir:$reference_dir $container \
    hisat2_extract_splice_sites.py $gtf_file > $reference_dir/Mus_musculus_splicesites.txt

apptainer exec --bind $reference_dir:$reference_dir $container \
    hisat2_extract_exons.py $gtf_file > $reference_dir/Mus_musculus_exons.txt


# providing the path where my reference genome for mus, gtf file, exon file and the splicesites file are located
# which is the same as reference_dir, so can use the same path or variable name  
#reference_genome_files="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome"

# Building HISAT2 index 
# --ss , it gives hisat2 a list of splice sites which we created earlier 
# --exon, it provides hisat2 with all exons that we craetedd earlier 
# -p is number of CPUs to be used for the job 
# hisat2-build , is the actual command that used inside the conatiner that will do the actual indexing  
#apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples $container \
#    hisat2-build \
#    --ss "$reference_dir/Mus_musculus_splicesites.txt" \
#    --exon "$reference_dir/Mus_musculus_exons.txt" \
#    -p ${SLURM_CPUS_PER_TASK} \
#    "$genome_fasta" "$index_prefix"


apptainer exec --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples $container \
    hisat2-build \
    -p ${SLURM_CPUS_PER_TASK} \
    "$genome_fasta" "$index_prefix"