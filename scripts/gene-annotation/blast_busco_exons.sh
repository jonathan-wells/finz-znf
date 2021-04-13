#!/usr/bin/env bash

for busco in 13830at7898; do
    tblastn \
        -query "../../data/danio-reads/${busco}.faa" \
        -subject ../../data/danio-reads/Danio_choprai_flanked_genes.fa \
        -outfmt "6 sseqid sstart send " \
        -evalue 1e-10 \
        -out test.out
    awk -v buscovar="${busco}" '{ OFS="\t" } {
        if ($2 > $3) 
            print $1, $3, $2, buscovar, ".",".", "BUSCO", "CDS", 0
        else
            print $1, $2, $3, buscovar, ".",".", "BUSCO", "CDS", 0
    }' test.out > test2.out 
done    
