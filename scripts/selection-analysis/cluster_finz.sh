#!/usr/bin/env bash

cd-hit \
    -i ../../data/seqs/cypriniformes_augustus_finz.fa \
    -o ../../data/cdhit_clusters \
    -d 100 \
    -M 0 \
    -T 12 \
    -sc 1 \
    -g 1 \
    -S 60 \
    -c 0.8
