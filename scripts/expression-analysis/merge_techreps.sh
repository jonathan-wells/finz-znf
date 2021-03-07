#!/usr/bin/env bash

while read line; do
    techrep=$(awk '{ print $1 }' <<< "$line")
    echo $techrep
    samtools sort \
        --threads 50 \
        "../../data/expression/STARaligned/${techrep}_Aligned.out.bam" \
        > "../../data/expression/STARaligned/${techrep}.sorted.bam"
done < ../../data/expression/White2017/accession_list.txt

while read line; do
    sample=$(awk '{ print $1 }' <<< "$line")
    techreps=$(awk '{ print $2, $3, $4, $5 }' <<< "$line")
    read -a array <<< "$techreps"
    echo $sample
    samtools merge \
        --threads 50 \
        -f \
        -h "../../data/expression/STARaligned/${array[0]}.sorted.bam" \
        -c \
        -p \
        "../../data/expression/STARaligned/${sample}.bam" \
        "../../data/expression/STARaligned/${array[0]}.sorted.bam" \
        "../../data/expression/STARaligned/${array[1]}.sorted.bam" \
        "../../data/expression/STARaligned/${array[2]}.sorted.bam" \
        "../../data/expression/STARaligned/${array[3]}.sorted.bam"
    samtools sort \
        -n \
        --threads 50 \
        "../../data/expression/STARaligned/${sample}.bam" \
        > "../../data/expression/STARaligned/${sample}.sorted.bam"
done < ../../data/expression/White2017/sample_list.txt
