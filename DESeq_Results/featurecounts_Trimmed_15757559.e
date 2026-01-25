
        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v2.0.6

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 15 BAM files                                     ||
||                                                                            ||
||                           SRR7821949_sorted.trimmed.bam                    ||
||                           SRR7821950_sorted.trimmed.bam                    ||
||                           SRR7821951_sorted.trimmed.bam                    ||
||                           SRR7821952_sorted.trimmed.bam                    ||
||                           SRR7821953_sorted.trimmed.bam                    ||
||                           SRR7821954_sorted.trimmed.bam                    ||
||                           SRR7821955_sorted.trimmed.bam                    ||
||                           SRR7821956_sorted.trimmed.bam                    ||
||                           SRR7821957_sorted.trimmed.bam                    ||
||                           SRR7821968_sorted.trimmed.bam                    ||
||                           SRR7821969_sorted.trimmed.bam                    ||
||                           SRR7821970_sorted.trimmed.bam                    ||
||                           SRR7821971_sorted.trimmed.bam                    ||
||                           SRR7821972_sorted.trimmed.bam                    ||
||                           SRR7821973_sorted.trimmed.bam                    ||
||                                                                            ||
||             Output file : gene_counts_trimmed.txt                          ||
||                 Summary : gene_counts_trimmed.txt.summary                  ||
||              Paired-end : yes                                              ||
||        Count read pairs : no                                               ||
||              Annotation : Mus_musculus.GRCm39.115.gtf (GTF)                ||
||      Dir for temp files : /data/users/dpatel/rna_seq_project_Bloodsamp ... ||
||                                                                            ||
||                 Threads : 4                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Mus_musculus.GRCm39.115.gtf ...                       ||
||    Features : 1299184                                                      ||
||    Meta-features : 78334                                                   ||
||    Chromosomes/contigs : 38                                                ||
||                                                                            ||
|| Process BAM file SRR7821949_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 66886669                                             ||
||    Successfully assigned alignments : 44102793 (65.9%)                     ||
||    Running time : 0.43 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821950_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 50255605                                             ||
||    Successfully assigned alignments : 33350646 (66.4%)                     ||
||    Running time : 0.25 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821951_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 74468684                                             ||
||    Successfully assigned alignments : 48903506 (65.7%)                     ||
||    Running time : 0.45 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821952_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 51921276                                             ||
||    Successfully assigned alignments : 34598052 (66.6%)                     ||
||    Running time : 0.26 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821953_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 48146575                                             ||
||    Successfully assigned alignments : 31678351 (65.8%)                     ||
||    Running time : 0.24 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821954_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 61005606                                             ||
||    Successfully assigned alignments : 33242130 (54.5%)                     ||
||    Running time : 0.27 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821955_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 79036596                                             ||
||    Successfully assigned alignments : 43579836 (55.1%)                     ||
||    Running time : 0.39 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821956_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 69330398                                             ||
||    Successfully assigned alignments : 35823257 (51.7%)                     ||
||    Running time : 0.36 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821957_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 62750937                                             ||
||    Successfully assigned alignments : 34074945 (54.3%)                     ||
||    Running time : 0.40 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821968_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 72394764                                             ||
||    Successfully assigned alignments : 39436758 (54.5%)                     ||
||    Running time : 0.34 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821969_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 42344068                                             ||
||    Successfully assigned alignments : 25243463 (59.6%)                     ||
||    Running time : 0.19 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821970_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 59592053                                             ||
||    Successfully assigned alignments : 32744739 (54.9%)                     ||
||    Running time : 0.26 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821971_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 58321369                                             ||
||    Successfully assigned alignments : 33322621 (57.1%)                     ||
||    Running time : 0.26 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821972_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 54950663                                             ||
||    Successfully assigned alignments : 33007718 (60.1%)                     ||
||    Running time : 0.24 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR7821973_sorted.trimmed.bam...                          ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 53597830                                             ||
||    Successfully assigned alignments : 34509935 (64.4%)                     ||
||    Running time : 0.24 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
|| Summary of counting results can be found in file "/data/users/dpatel/rna_  ||
|| seq_project_Bloodsamples/Trimmed_data/countsdata_trimmed/gene_counts_trim  ||
|| med.txt.summary"                                                           ||
||                                                                            ||
\\============================================================================//

