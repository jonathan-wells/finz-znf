#!/usr/bin/env bash

# First activate busco venv with:
# $ conda activate buscoenv

genomedir="/workdir/jnw72/Genomes/ncbi-genomes-2020-07-09"
declare -A genomes
# while read line; do
#     key=`awk '{ print $1 }' <<< $line`
#     data=`awk '{ print $2 }' <<< $line`
#     genomes[$key]="$data"
# done < ../../data/species_genomes.txt

genomes=(
    ['Labeo_gonius']='GCA_013461565.1_ASM1346156v1_genomic.fna'
    ['Paedocypris_carbunculus']='QDDN_carbunculus.CA_carbunculus.softmasked.fna'
    ['Paedocypris_micromegethes']='QDDN_micromegethes.CA_micromegethes.softmasked.fna'
    )


for species in ${!genomes[@]}; do
    busco \
        -i "${genomedir}/${genomes[$species]}" \
        -l "actinopterygii_odb10" \
        -o $species \
        --out_path "../../data/busco-out/" \
        -m "genome" \
        -c 40 \
        -f
done
