#!/usr/bin/env bash

# First activate busco venv with:
# $ conda activate buscoenv

genomedir="/workdir/jnw72/Genomes"
declare -A genomes
# while read line; do
#     key=`awk '{ print $1 }' <<< $line`
#     data=`awk '{ print $2 }' <<< $line`
#     genomes[$key]="$data"
# done < ../../data/species_genomes.txt

genomes=(
[Puntigrus_tetrazona]=GCA_018831695.1_ASM1883169v1_genomic.fna
[Sinocyclocheilus_maitianheensis]=GCA_018148995.1_ASM1814899v1_genomic.fna
[Sinocyclocheilus_anophthalmus]=GCA_018155175.1_ASM1815517v1_genomic.fna
[Gobiocypris_rarus]=GCA_018491645.1_ASM1849164v1_genomic.fna
[Paracanthobrama_guichenoti]=GCA_018749465.1_ASM1874946v1_genomic.fna
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
