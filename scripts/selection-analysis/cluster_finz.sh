#!/usr/bin/env bash

for id in 75 80 85 90 95; do
    cd-hit \
        -i ../../data/seqs/cypriniformes_augustus_finz.fa \
        -o "../../data/cdhit-out/finz_percentid_${id}" \
        -d 100 \
        -M 0 \
        -T 12 \
        -sc 1 \
        -g 1 \
        -S 60 \
        -c "0.${id}"

    for i in {0..15}; do
        ./parse_cdhit_clusters.py getseqs \
            -f "../../data/cdhit-out/finz_percentid_${id}.clstr" \
            -i $i \
            -c True > "../../data/cdhit-out/finz_percentid_${id}_c${i}.fna" 
    done

done

