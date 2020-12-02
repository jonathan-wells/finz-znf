#!/usr/bin/env bash

while read line; do
    sample=$(awk '{ print $1 }' <<< "$line")
    echo $sample
    TEcount \
        -b "../../data/expression/STARaligned/${sample}.bam" \
        --GTF ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf \
        --TE ../../data/gffs/danRer11.nonalt.tetranscripts.gtf \
        --sortByPos \
        --project "../../data/expression/TEcount-out/${sample}" &
done < ../../data/expression/White2017/sample_list.txt

