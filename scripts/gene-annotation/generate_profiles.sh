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

# For genes used to generate species phylogeny
# These alignments were downloaded from EggNog on 
for prot in FICD GLYT ILF2 MYH6 RAG1 RAG2 TBR1; do
    echo "Generating $prot profiles"
    
    hmmbuild "../../data/phmms/${prot}.hmm" "../../data/phmms/${prot}.fa"
    hmmemit -c "../../data/phmms/${prot}.hmm" > "../../data/phmms/${prot}.consensus.fa"

    prepareAlign \
        < "../../data/phmms/${prot}.fa" \
        > "../../data/phmms/${prot}.prepped.fa"
    
    /usr/local/Cellar/augustus/3.3.3_1/scripts/msa2prfl.pl \
        --qij=/usr/local/Cellar/augustus/3.3.3_1/config/profile/blosum62.qij \
        "../../data/phmms/${prot}.prepped.fa" \
        > "../../data/phmms/${prot}.prfl"
done
