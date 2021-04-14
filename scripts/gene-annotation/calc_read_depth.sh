#!/usr/bin/env bash

declare -A species
species=(
    [Danio_albolineatus]=ERR3284960
    [Danio_jaintianensis]=ERR3284968
    [Danio_choprai]=ERR3284956
    [Danio_tinwini]=ERR3284972
    )

for sp in ${!species[@]}; do
    samtools sort \
        -O BAM \
        -o "../../data/danio-reads/${sp}_${species[$sp]}.sorted.bam" \
        -@ 30 \
        "../../data/danio-reads/${sp}_${species[$sp]}.sam"
    
    samtools depth \
        -o "../../data/danio-reads/${sp}_depth.out" \
        --reference "../../data/danio-reads/${sp}_flanked_genes.fa" \
        -b "../../data/danio-reads/${sp}_busco_cds.bed" \
        "../../data/danio-reads/${sp}_${species[$sp]}.sorted.bam"

done
