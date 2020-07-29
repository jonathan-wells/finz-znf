#!/usr/bin/env bash

################################################################################
## Use BLAST to search for finz domain in translated genome sequences 
################################################################################

for file in /Users/jonwells/Genomes/ncbi-genomes-2020-07-09/*.fna; do
    if [[ $file =~ /Users/jonwells/Genomes/ncbi-genomes-2020-07-09/(GC[A-Z]_[0-9]+\.[0-9]) ]]; then
        acc=${BASH_REMATCH[1]}
        spec=`datasets assembly-descriptors accession $acc | 
            python -c 'import sys, json; print(json.load(sys.stdin)["datasets"][0]["org"]["sci_name"].replace(" ", "_"))'`
        echo $spec
        tblastn \
            -query ../data/seqs/finz_domain_consensus.fa \
            -subject $file \
            -outfmt 6 \
            -out "../data/blast-out/${spec}_finz.out"

        bioawk -c fastx '{ 
            OFS="\t" 
        } { 
            print $name, length($seq) 
        }' $file > "../data/beds/${spec}.genome" 
        
    fi
done
