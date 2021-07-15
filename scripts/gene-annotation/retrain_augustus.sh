#!/usr/bin/env bash

/usr/local/Augustus/scripts/gff2gbSmallDNA.pl \
    ../../data/gffs/ensembl_training_finz.gff \
    /Users/jonwells/Genomes/Cypriniformes/GCF_000002035.6_GRCz11_genomic.nonalt.fna \
    6000 \
    ../../data/gffs/ensembl_training_finz.gb

/usr/local/Augustus/scripts/randomSplit.pl \
    ../../data/gffs/ensembl_training_finz.gb 20

/usr/local/Augustus/bin/etraining \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.train

/usr/local/Augustus/bin/augustus \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.test \
    > initial_eval.out \
    2> initial_eval.err

/usr/local/Augustus/scripts/optimize_augustus.pl \
    --species=zebrafish_new \
    --cpus=8 \
    --rounds=1 \
    --UTR=on \
    ../../data/gffs/ensembl_finz_znf_expressed.gb.train

/usr/local/Augustus/bin/augustus \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.test \
    > final_eval.out \
    2> final_eval.err

