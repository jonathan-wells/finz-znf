#!/usr/bin/env bash

for prot in CDAN1 DCST2 FICD GLYT ILF2 RAG1 RAG2 RWDD3; do
    ginsi "../../data/species-phylogeny/multispecies_${prot}.fa" \
        > "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" 
    trimal \
        -automated1 \
        -in "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" \
        -out "../../data/species-phylogeny/multispecies_${prot}.trimmed.fa"
done

for prot in TPO FSHR WFIKKN1; do
    ginsi "../../data/species-phylogeny/of_multispecies_${prot}.fa" \
        > "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" 
    trimal \
        -automated1 \
        -in "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" \
        -out "../../data/species-phylogeny/multispecies_${prot}.trimmed.fa"
done

