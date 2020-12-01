#!/usr/bin/env bash

# Generate genome indices
/programs/STAR-2.7.5a/bin/Linux_x86_64/STAR \
    --genomeDir "../../data/expression/STARgenome" \
    --genomeFastaFiles "../../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna" \
    --runMode genomeGenerate \
    --runThreadN 25 \
    --sjdbGTFfile "../../data/expression/STARgenome/danRer11.nonalt.genes_tes.gtf" \
    --sjdbOverhang 99 

# Align reads
for readfile in $(cat ../../data/expression/White2017/accession_list.txt); do
    echo Processing $readfile ...
    /programs/STAR-2.7.5a/bin/Linux_x86_64/STAR \
        --genomeDir "../../data/expression/STARgenome" \
        --readFilesIn "../../data/expression/White2017/${readfile}_1.fastq.gz" "../../data/expression/White2017/${readfile}_2.fastq.gz" \
        --readFilesCommand gunzip -c \
        --runThreadN 8 \
        --outSAMtype BAM Unsorted \
        --runMode alignReads \
        --outFilterMultimapNmax 5000 \
        --winAnchorMultimapNmax 5000 \
        --outMultimapperOrder Random \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix "../../data/expression/STARaligned/${readfile}_"
done
