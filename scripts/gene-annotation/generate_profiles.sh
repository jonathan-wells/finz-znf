#!/usr/bin/env bash

# For finz-znfs specifically
echo "Generating FINZ-ZNF profile"
prepareAlign \
    < ../../data/phmms/drerio_finz_expressed.best.fas \
    > ../../data/phmms/drerio_finz_expressed.prepped.fa

/usr/local/Cellar/augustus/3.3.3_1/scripts/msa2prfl.pl \
    --qij=/usr/local/Cellar/augustus/3.3.3_1/config/profile/blosum62.qij \
    ../../data/phmms/drerio_finz_expressed.prepped.fa \
    > ../../data/phmms/drerio_finz_expressed.prfl

# NOTE: The above commands were used on a set of robustly expressed D. rerio
# FiNZ-ZnF genes previously annotated by Ensembl and White et al. (2017). Since
# this is inherently biased towards D. rerio, a second iteration of gene
# annotation was carried out using a profile generated from multiple sequence
# alignments of cypriniforme ZnFs and FiNZ domains from a larger subset of
# species. These profiles (available in data dir) were used for final analyses.
