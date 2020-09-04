#!/usr/bin/env bash

for prot in FICD GLYT ILF2 RAG1 RAG2; do
    ginsi "../../data/species-phylogeny/multispecies_${prot}.fa" \
        > "../../data/species-phylogeny/multispecies_${prot}.ginsi.fa" 
done


