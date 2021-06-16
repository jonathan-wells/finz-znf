#!/usr/bin/env bash

# # Extract biological replicates corresponding to stages where FiNZ-ZnF
# # expression peaks.
# awk '{
#     if ($NF == "RNASeq" && 
#         ($6 == "Dome" || 
#          $6 == "128-cell" || 
#          $6 == "1k-cell" || 
#          $6 == "50pc-epiboly" || 
#          $6 == "Shield" || 
#          $6 == "75pc-epiboly"))
#     print "../../data/expression/STARaligned/"$3".bam";
# }' ../../data/expression/White2017/elife-30860-supp1-v1.tsv > merge_bams.names

# # Merge and sort relevant bamfiles
# bamtools merge \
#     -list merge_bams.names \
#     -out peak_stages.bam
# samtools sort \
#     -o peak_stages.sorted.bam \
#     -O BAM \
#     --reference ../../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna \
#     -@ 25 \
#     peak_stages.bam

# Run Trinity in genome-guided mode to extract transcripts
/workdir/jnw72/Software/trinityrnaseq-v2.12.0/Trinity \
    --genome_guided_bam peak_stages.sorted.bam \
    --max_memory 250G \
    --genome_guided_max_intron 10000 \
    --CPU 20 \
    --output ../../data/expression/trinity-transcripts
