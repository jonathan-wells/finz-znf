#!/usr/bin/env bash

genomedir="/local/workdir/jnw72/Genomes/ncbi-genomes-2020-07-09"
declare -A genomes

while read line; do
    key=`awk '{ print $1 }' <<< $line`
    data=`awk '{ print $2 }' <<< $line`
    genomes[$key]="$data"
done < ../../data/species_genomes.txt

buscodir="../../data/busco-out"
buscooutdir="../../data/busco-out/super"

rm "${buscooutdir}/*"

scseqs="run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences"
while read -r buscoid; do
    for species in ${!genomes[@]}; do
        buscoseq="${buscodir}/${species}/${scseqs}/${buscoid}.faa"
        sed "s/>/>${species}_${buscoid}_/" $buscoseq > tmp.fa
        cat tmp.fa >> "${buscooutdir}/${buscoid}.fa"
    done
done < ../../data/busco_ids.txt
rm tmp.fa