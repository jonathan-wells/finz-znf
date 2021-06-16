#!/usr/bin/env bash

while read line; do
    sample=$(awk '{ print $1 }' <<< "$line")
    echo $sample
    samtools view -h "../../data/expression/STARaligned/${sample}.sorted.bam" |
        awk 'substr($0,1,1)=="@" || ($9<=5000 && $9>=-5000)' |
        samtools view -b > "../../data/expression/STARaligned/${sample}.inslt5kb.bam" &
    sleep 5
done < ../../data/expression/White2017/sample_list.txt
