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
    ['Labeo_gonius']='GCA_014525385.1_HRRL_MCC_1.0_genomic.fna'
    ['Chanos_chanos']='GCF_902362185.1_fChaCha1.1_genomic.fna'
    ['Triplophysa_dalaica']='GCA_015846415.1_ASM1584641v1_genomic.fna'
)

for species in ${!genomes[@]}; do
    busco \
        -i "${genomedir}/${genomes[$species]}" \
        -l "actinopterygii_odb10" \
        -o $species \
        --out_path "../../data/busco-out/" \
        -m "genome" \
        -c 40
done
