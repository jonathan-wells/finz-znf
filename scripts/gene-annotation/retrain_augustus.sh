#!/usr/bin/env bash

/usr/local/Cellar/augustus/3.3.3_2/scripts/gff2gbSmallDNA.pl \
    ../../data/gffs/ensembl_training_finz.gff \
    /Users/jonwells/Genomes/Cypriniformes/GCF_000002035.6_GRCz11_genomic.nonalt.fna \
    6000 \
    ../../data/gffs/ensembl_training_finz.gb

/usr/local/Cellar/augustus/3.3.3_2/scripts/randomSplit.pl \
    ../../data/gffs/ensembl_training_finz.gb 12

/usr/local/Cellar/augustus/3.3.3_2/bin/etraining \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.train

/usr/local/Cellar/augustus/3.3.3_2/bin/augustus \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.test \
    > initial_eval.out \
    2> initial_eval.err

# /usr/local/Cellar/augustus/3.3.3_2/scripts/optimize_augustus.pl \
#     --species=zebrafish_new \
#     --cpus=8 \
#     --rounds=4 \
#     --UTR=on \
#     ../../data/gffs/ensembl_training_finz.gb.train

# /usr/local/Cellar/augustus/3.3.3_2/scripts/optimize_augustus.pl \
#     --species=zebrafish_new \
#     --cpus=8 \
#     --rounds=3 \
#     --trainOnlyUtr=1 \
#     --UTR=on \
#     ../../data/gffs/ensembl_training_finz.gb.train

/usr/local/Cellar/augustus/3.3.3_2/bin/augustus \
    --species=zebrafish_new \
    --UTR=on \
    ../../data/gffs/ensembl_training_finz.gb.test \
    > final_eval.out \
    2> final_eval.err

