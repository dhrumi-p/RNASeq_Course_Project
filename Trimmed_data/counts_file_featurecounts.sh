#!/bin/bash

#SBATCH --job-name=featureCounts_Trimmed
#SBATCH --cpus-per-task=6
#SBATCH --time=05:00:00
#SBATCH --mem=16GB
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/output_trimmed/featurecounts_Trimmed_%j.o
#SBATCH --error=/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/error_trimmed/featurecounts_Trimmed_%j.e
#SBATCH --mail-user=dhrumi.patel@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

# Paths
INPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/mapped_files_trimmed"
OUTPUT_DIR="/data/users/dpatel/rna_seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed"
GTF_FILE="/data/users/dpatel/rna_seq_project_Bloodsamples/Mouse_reference_genome/Mus_musculus.GRCm39.115.gtf"
CONTAINER="/containers/apptainer/subread_2.0.6.sif" # feature counts 2.0.6 (version)

# ADD THIS LINE - Define the output file
COUNTS_OUTPUT="$OUTPUT_DIR/gene_counts_trimmed.txt"

# feature counts will only save things to output file when there is some file mentioned to ti and thus teh results
# wont go to directory as such 

# Collect BAM files
BAM_FILES=()
for f in "$INPUT_DIR"/*_sorted.trimmed.bam; do
    [[ -f "$f" ]] && BAM_FILES+=("$f")
done

apptainer exec \
    --bind /data/users/dpatel/rna_seq_project_Bloodsamples:/data/users/dpatel/rna_seq_project_Bloodsamples \
    "$CONTAINER" \
    featureCounts \
        -T 4 \
        -s 2 \
        -p \
        -a "$GTF_FILE" \
        -o "$COUNTS_OUTPUT" \
        -t exon \
        -g gene_id \
        "${BAM_FILES[@]}"
   
#   -T 4 \  # cpus-per-task
# -s is the strandedness of the RNA seq library, default ia 0 meaning unstranded
# -s 1 meaning reads are from sense starnd and s -2 means reads are from antisense strand
#   -p \ # it takes in to account pair end reads 
#    -a "$GTF_FILE" \  # Annotation file which contains all the information about the exons,introns and other stuff in the mouse reference genome 
#    -o "$COUNTS_OUTPUT" \  # output file where the counts tabel will be stored 
#    -t exon \   # -t tells featurecounts which feature to take in consideration to match reads, here it will be exon
#    -g gene_id \  # it is telling that gene_id to be used from GTF file, to group the reads 
#    "${BAM_FILES[@]}" # take each sorted.bam.bai files and pass it as argument to featurecounts individually 



# -t shows what exactly we are tring to map/calculate or in other words against what the reads needs to be mapped 
# t gave the information on where to map the reads 
# -g will give information how to club the reads/ taking what to consideration and give one number (basically on what basis to group the reads with)
# since here we would want how many reads aligned to each genes coding for RNA(transcript) (exon part of genes)
# and thus we have taken exons in -t and gene_id in -g 
# suppose, if we wish to group reads using gene_name, then in -g instead of gene_is, we can use gene_name