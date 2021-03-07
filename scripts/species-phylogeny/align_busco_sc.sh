#!/usr/bin/env bash

buscodir="../../data/busco-out/super"
while read -r buscoid; do
    ginsi "${buscodir}/${buscoid}.fa" \
        > "../../data/species-phylogeny/multispecies_${buscoid}.ginsi.fa" 
    trimal \
        -automated1 \
        -in "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" \
        -out "../../data/species-phylogeny/multispecies_${prot}.trimmed.fa"
done < ../../data/busco_ids.txt

