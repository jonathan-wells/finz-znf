#!/usr/bin/env bash

species=(
    Danio_choprai
    Danio_jaintianensis
    Danio_albolineatus
    Danio_tinwini
    )

for sp in ${species[@]}; do
    echo $sp
    rg "\d+at\d+" ../../data/danio-reads/${sp}_flanked_genes.fa |
        sed "s/>${sp}_//" > tmp_buscos.txt

    while read -r busco; do
        tblastn \
            -query "../../data/busco-out/${sp}/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/${busco}.faa" \
            -subject ../../data/danio-reads/${sp}_flanked_genes.fa \
            -outfmt "6 sseqid sstart send " \
            -evalue 1e-10 \
            -out tmp.out
        awk -v buscovar="${busco}" '{ OFS="\t" } {
            if ($2 > $3) 
                print $1, $3, $2, buscovar, ".",".", "BUSCO", "CDS", 0
            else
                print $1, $2, $3, buscovar, ".",".", "BUSCO", "CDS", 0
        }' tmp.out >> "${sp}.out" 
    done < tmp_buscos.txt 
done 
rm tmp_buscos.txt tmp.out

