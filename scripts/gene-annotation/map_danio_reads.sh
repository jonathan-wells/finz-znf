#!/usr/bin/env bash

declare -A readaccs
readaccs=(
    [Danio_albolineatus]=ERR3284960
    [Danio_jaintianensis]=ERR3284968
    [Danio_choprai]=ERR3284956
    [Danio_tinwini]=ERR3284972
)

genomedir="/local/workdir/Genomes/Cypriniformes/ncbi-genomes-2020-07-09"
declare -A genomes
genomes=(
    [Danio_albolineatus]="${genomedir}/GCA_903798035.1_fDanAlb1.1_genomic.fna"
    [Danio_jaintianensis]="${genomedir}/GCA_903798115.1_fDanJai1.1_genomic.fna"
    [Danio_choprai]="${genomedir}/GCA_903798125.1_fDanCho1.1_genomic.fna"
    [Danio_tinwini]="${genomedir}/GCA_903798205.1_fDanTin1.1_genomic.fna"
    )

for species in ${!genomes[@]}; do
    /programs/minimap2-2.17/minimap2 \
        -ax sr \
        -t 20 \
        --secondary=no \
        ../../data/danio-reads/${species}_flanked_genes.fa \
        ../../data/danio-reads/${species}_${readaccs[$species]}_1.fastq \
        ../../data/danio-reads/${species}_${readaccs[$species]}_2.fastq \
        > ../../data/danio-reads/${species}_${readaccs[$species]}.sam
done
