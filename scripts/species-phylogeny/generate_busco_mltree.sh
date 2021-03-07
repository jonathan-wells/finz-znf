#!/usr/bin/env bash

treedir="../../data/species-phylogeny"
iqtree2 \
    -s "${treedir}/busco_supermatrix.fa" \
    -p "${treedir}/busco_supermatrix_partitions.nex" \
    -T 40 \
    -B 5000
