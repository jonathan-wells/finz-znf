#!/usr/bin/env bash

/usr/local/Cellar/augustus/3.3.3_1/scripts/autoAug.pl \
    -g ~/Desktop/Chr4.fa \
    -t ../../data/gffs/ensembl_finz_znf.gff \
    --species=DanioChr4 \
    -v 5 \
    --noninteractive


