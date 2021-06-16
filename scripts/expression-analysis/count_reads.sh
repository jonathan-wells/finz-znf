#!/usr/bin/env bash

while read line; do
    sample=$(awk '{ print $1 }' <<< "$line")
    echo $sample
    /programs/TEtranscripts-2.2.1/bin/TEcount \
        -b "../../data/expression/STARaligned/${sample}.inslt5kb.bam" \
        --GTF ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf \
        --TE ../../data/gffs/danRer11.nonalt.tetranscripts.gtf \
        --mode multi \
        --format BAM \
        --stranded reverse \
        --project "../../data/expression/TEcount-out/${sample}.inslt5kb" &
    sleep 30 
done < ../../data/expression/White2017/sample_list.txt

