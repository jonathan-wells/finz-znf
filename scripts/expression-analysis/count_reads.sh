#!/usr/bin/env bash

for readfile in $(cat ../../data/expression/White2017/accession_list.txt); do
    echo Processing $readfile ...
    TEcount \
        -b "../../data/expression/STARaligned/${readfile}_Aligned.out.bam" \
        --GTF ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf \
        --TE ../../data/gffs/danRer11.nonalt.tetranscripts.gtf \
        --project "../../data/expression/TEcount-out/${readfile}"
done
