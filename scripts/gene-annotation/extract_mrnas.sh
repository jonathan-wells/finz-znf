#!/usr/bin/env bash

# Extract biological replicates corresponding to stages where FiNZ-ZnF
# expression peaks.

for stage in "Dome" "128-cell" "1k-cell" "50pc-epiboly" "Shield" "75pc-epiboly"; do
    echo extracting $stage stage
    awk -v stage=$stage '{
        if ($NF == "RNASeq" && $6 == stage ) 
            print "../../data/expression/STARaligned/"$3".bam";
    }' ../../data/expression/White2017/elife-30860-supp1-v1.tsv > "${stage}.names"

    # Merge and sort relevant bamfiles
    samtools merge \
        -b "${stage}.names" "${stage}.bam" \
        -h \
        -c \
        -p \
        -@ 25
    samtools sort \
        -o "${stage}.sorted.bam" \
        -O BAM \
        --reference ../../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna \
        -@ 25 \
        "${stage}.bam"
    rm "${stage}.bam" "${stage}.names"


    # Run Trinity in genome-guided mode to extract transcripts
    echo running trinity for stage $stage
    /workdir/jnw72/Software/trinityrnaseq-v2.12.0/Trinity \
        --genome_guided_bam "${stage}.sorted.bam"  \
        --max_memory 500G \
        --genome_guided_max_intron 10000 \
        --CPU 25 \
        --output "../../data/expression/trinity-transcripts/trinity-${stage}"
done
