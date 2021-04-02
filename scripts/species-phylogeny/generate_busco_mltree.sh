#!/usr/bin/env bash

treedir="../../data/species-phylogeny/iqtree-data"
iqtree2 \
    -s "${treedir}/busco_supermatrix.fa" \
    -p "${treedir}/busco_supermatrix_partitions.nex" \
    --date "${treedir}/cypriniformes_datefile.txt" \
    --date-tip 0 \
    --date-ci 100 \
    --sampling "GENESITE" \
    --boot-trees \
    --alrt 1000 \
    -B 5000 \
    -T 40 
