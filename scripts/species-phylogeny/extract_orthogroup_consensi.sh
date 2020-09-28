#!/usr/bin/env bash

orthodir="../../data/species-phylogeny"
for file in ${orthodir}/Results_Sep08/Single_Copy_Orthologue_Sequences/*; do
    og="${file##*/}"
    mafft $file > "${orthodir}/orthogroup-consensi/${og}"
    hmmbuild tmp.hmm "${orthodir}/orthogroup-consensi/${og}"
    hmmemit -c tmp.hmm > "${orthodir}/orthogroup-consensi/${og}"
done

rm tmp.hmm
