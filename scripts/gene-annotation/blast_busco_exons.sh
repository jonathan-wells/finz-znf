#!/usr/bin/env bash

species=(
    Danio_choprai
    Danio_jaintianensis
    Danio_albolineatus
    Danio_tinwini
    )

for sp in ${species[@]}; do
    echo $sp
    rm -r "tmp_${sp}"
    mkdir "tmp_${sp}"
    rg "\d+at\d+" "../../data/danio-reads/${sp}_flanked_genes.fa" |
        cut -c 2- > "tmp_${sp}/tmp_buscos.txt"
    seqtk subseq "../../data/danio-reads/${sp}_flanked_genes.fa" "tmp_${sp}/tmp_buscos.txt" \
        > "tmp_${sp}/busco_seqs.fa"
    fastaexplode --fasta "tmp_${sp}/busco_seqs.fa" -d "tmp_${sp}"
    
    # Strip species names
    sed "s/${sp}_//" "tmp_${sp}/tmp_buscos.txt" > tmp.txt
    mv tmp.txt "tmp_${sp}/tmp_buscos.txt"
    
    echo "blasting seqs"
    while read -r busco; do
        tblastn \
            -query "../../data/busco-out/${sp}/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/${busco}.faa" \
            -subject "tmp_${sp}/${sp}_${busco}.fa" \
            -outfmt "6 sseqid sstart send " \
            -evalue 1e-10 \
            -out tmp.out
        awk -v buscovar="${busco}" '{ OFS="\t" } {
            if ($2 > $3) 
                print $1, $3, $2, buscovar, ".",".", "BUSCO", "CDS", 0
            else
                print $1, $2, $3, buscovar, ".",".", "BUSCO", "CDS", 0
        }' tmp.out >> "${sp}.out" 
    done < "tmp_${sp}/tmp_buscos.txt"
done 
rm -r tmp*

