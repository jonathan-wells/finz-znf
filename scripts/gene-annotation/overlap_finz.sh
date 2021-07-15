#!/usr/bin/env bash
 
# This script finds canonical finz proteins in the RefSeq and Ensembl gene sets, 
# based on shared presence of finz and zinc finger domains. It subsets these
# genes from the corresponding gff files in order to calculate degree of overlap
# the new predictions from augustus-ppx.

# TODO: change denovo_finz block to use existing Danio rerio annotations.

declare -A annotations
annotations=(
    [RefSeq]=../../data/seqs/GCF_000002035.6_GRCz11_protein.faa
    [Ensembl]=../../data/seqs/Danio_rerio.GRCz11.pep.all.fa
)

for i in ${!annotations[@]}; do
    echo Searching $i
    hmmsearch \
        --tblout tmp.out \
        -E 1e-04 \
        ../../data/phmms/finz_seed.hmm \
        ${annotations[$i]}
    rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names
    
    hmmsearch \
        --tblout tmp.out \
        -E 1e-04 \
        ../../data/phmms/PF00096_seed.hmm \
        ${annotations[$i]}
    rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

    cat finz.names c2h2.names | sort | uniq -d > "${i}_finz_znf.names"
    ./extract_gene_names.py $i > "../../data/gffs/${i}_finz_znf.gff"

done
./rename_ensembl_chroms.py ../../data/gffs/ensembl_finz_znf.gff ../../data/gffs/ensembl_finz_znf.gff

rg Danio_rerio ../../data/seqs/cypriniformes_augustus_finz.fa |
    cut -c 14- > denovo_finz_znf.names
./extract_gene_names.py denovo > "../../data/gffs/denovo_finz_znf.gff"
rm *.names tmp.out

# Then calculate overlap
echo calculating overlap
rg "\tmRNA\t" ../../data/gffs/refseq_finz_znf.gff \
    | gff2bed > ../../data/beds/refseq_finz_znf.transcripts.bed
rg "\tmRNA\t" ../../data/gffs/ensembl_finz_znf.gff \
    | gff2bed > ../../data/beds/ensembl_finz_znf.transcripts.bed
rg "\ttranscript\t" ../../data/gffs/denovo_finz_znf.gff \
    | gff2bed > ../../data/beds/denovo_finz_znf.transcripts.bed

bedtools intersect \
    -a ../../data/beds/Ensembl_finz_znf.transcripts.bed \
    -b ../../data/beds/RefSeq_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.7 \
    -r \
    > ensembl_refseq.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ../../data/beds/RefSeq_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.7 \
    -r \
    > denovo_refseq.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ../../data/beds/Ensembl_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.7 \
    -r \
    > denovo_ensembl.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ensembl_refseq.bed \
    -wa \
    -wb \
    -s \
    -f 0.7 \
    -r \
    > denovo_ensembl_refseq.bed

./count_genes.py > "../../data/finz_znf_overlap_70.txt"
rm *.bed
