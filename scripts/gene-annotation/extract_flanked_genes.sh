#!/usr/bin/env bash

# This script extracts fasta sequence from around the gene bodies of predicted
# finz-znfs and a set of 10 busco single-copy orthologues. These sequences are
# then used as a reference on which to align raw genomic reads. The purpose of
# this is to test the completeness of the set of predicted finz genes. If
# coverage is substantially higher over the finz sequences than over the
# single-copy busco genes, that could indicate that finz genes have been missed
# by augustus.
# 

declare -A readaccs
readaccs=(
    [Danio_albolineatus]=ERR3284960
    [Danio_jaintianensis]=ERR3284968
    [Danio_choprai]=ERR3284956
    [Danio_tinwini]=ERR3284972
)

genomedir="/Users/jonwells/Genomes/Cypriniformes/ncbi-genomes-2020-07-09"
declare -A genomes
genomes=(
    [Danio_albolineatus]="${genomedir}/GCA_903798035.1_fDanAlb1.1_genomic.fna"
    [Danio_jaintianensis]="${genomedir}/GCA_903798115.1_fDanJai1.1_genomic.fna"
    [Danio_choprai]="${genomedir}/GCA_903798125.1_fDanCho1.1_genomic.fna"
    [Danio_tinwini]="${genomedir}/GCA_903798205.1_fDanTin1.1_genomic.fna"
    )

seqfile="../../data/seqs/cypriniformes_augustus_finz.fa"
for species in ${!readaccs[@]}; do
    echo $species
    
    # Generate genome coords file
    bioawk -c fastx '{ OFS="\t"} { print $name, length($seq) }' \
        "../../data/seqs/${species}_finz_blocks.fa" > blocks_genome.bed
    
    # Extract finz genes from gff
    rg $species $seqfile |
        sed "s/>${species}_/ID=/" |
        sed 's/\.t1/$/' > genepattern.txt
    rg -f genepattern.txt "../../data/gffs/${species}_augustus_finz.gff" > tmp.gff
    
    sed 's/\$/($|\\\.)/' genepattern.txt > genepattern2.txt
    rg -f genepattern2.txt "../../data/gffs/${species}_augustus_finz.gff" > tmp2.gff
    ./offset_gffs.py tmp2.gff tmp3.gff 200 "${species}_"  
    gff2bed < tmp3.gff |
        rg "\tCDS\t" > "../../data/danio-reads/${species}_busco_cds.bed" 
    
    # Extend flanking regions to each gene
    bedtools slop \
        -i tmp.gff \
        -g blocks_genome.bed \
        -b 200 | 
        awk '{ OFS="\t" } { print $1, $2, $9, $4, $5, $6, $7, $8 }' |
        sed "s/ID=/${species}_/" > tmp.bed

    # Extract finz fasta seqs
    bedtools getfasta \
        -fi "../../data/seqs/${species}_finz_blocks.fa" \
        -bed tmp.bed \
        -fo "../../data/danio-reads/${species}_flanked_genes.fa" \
        -nameOnly

    # Extract busco fasta seqs
    bedtools getfasta \
        -fi ${genomes[$species]} \
        -bed "../../data/danio-reads/${species}_busco.bed" \
        -nameOnly >> "../../data/danio-reads/${species}_flanked_genes.fa" 
done

rm tmp.bed tmp{,2,3}.gff genepattern.txt genepattern2.txt blocks_genome.bed
