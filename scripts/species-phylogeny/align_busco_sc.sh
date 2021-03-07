#!/usr/bin/env bash

buscodir="../../data/busco-out/super"
outdir="../../data/species-phylogeny/aligned-busco"
while read -r buscoid; do
    /programs/mafft/bin/ginsi \
        --thread 20 \
        "${buscodir}/${buscoid}.fa" > "${outdir}/${buscoid}.ginsi.fa" 
    /programs/trimal-1.4/source/trimal \
        -automated1 \
        -in "${outdir}/${buscoid}.ginsi.fa" \
        -out "${outdir}/${buscoid}.trimmed.fa"
done < ../../data/busco_ids.txt

