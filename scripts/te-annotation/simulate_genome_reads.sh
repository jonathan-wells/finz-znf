#!/usr/bin/env bash

genomedir="/local/workdir/Genomes/Cypriniformes"
declare -A genomes
while read line; do
    key=`awk '{ print $1 }' <<< $line`
    data=`awk '{ print $2 }' <<< $line`
    genomes[$key]="$data"
done < ../../data/species_genomes.txt

for species in ${!genomes[@]}; do
    genomefile=${genomes[$species]}
    /programs/art_bin_MountRainier/art_illumina \
        -ss HS25 \
        -ef \
        -d "${species}-"\
        -i "${genomedir}/${genomefile}" \
        -l 150 \
        -f 10 \
        -na \
        -o "../../data/simulated-reads/${species}" &

done
