#!/usr/bin/env bash

./njtree_labels.py ../../data/seqs/cypriniformes_augustus_finz.fa tmp.fas

# java \
#     -Xmx16g \
#     -jar /workdir/jnw72/Software/jD2Stat_1.0_SourceCode/jD2Stat_1.0.jar \
#     -d all \
#     -t 20 \
#     -a aa \
#     -i tmp.fas \
#     -o ../../data/selection-analysis/cypriniformes_augustus_finz_D2k8.mat

java \
    -Xmx60g \
    -jar /workdir/jnw72/Software/jD2Stat_1.0_SourceCode/jD2Stat_1.0.jar \
    -t 60 \
    -a aa \
    -n 1 \
    -i tmp.fas \
    -o ../../data/selection-analysis/cypriniformes_augustus_finz_D2nk8.mat

java \
    -Xmx60g \
    -jar /workdir/jnw72/Software/jD2Stat_1.0_SourceCode/jD2Stat_1.0.jar \
    -d all \
    -t 60 \
    -k 4 \
    -a aa \
    -i tmp.fas \
    -o ../../data/selection-analysis/cypriniformes_augustus_finz_D2k4.mat

# java \
#     -Xmx60g \
#     -jar /workdir/jnw72/Software/jD2Stat_1.0_SourceCode/jD2Stat_1.0.jar \
#     -t 40 \
#     -a aa \
#     -k 4 \
#     -n 1 \
#     -i tmp.fas \
#     -o ../../data/selection-analysis/cypriniformes_augustus_finz_D2nk4.mat
rm tmp.fas
