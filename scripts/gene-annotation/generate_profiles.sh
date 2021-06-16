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

