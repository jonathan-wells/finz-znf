#!/usr/bin/env bash
 
# This script finds canonical finz proteins in the RefSeq and Ensembl gene sets, 
# based on shared presence of finz and zinc finger domains. It subsets these
# genes from the corresponding gff files in order to calculate degree of overlap
# the new predictions from augustus-ppx.

# For RefSeq
hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/finz_seed.hmm \
    ../../data/seqs/GCF_000002035.6_GRCz11_protein.faa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/PF00096_seed.hmm \
    ../../data/seqs/GCF_000002035.6_GRCz11_protein.faa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

cat finz.names c2h2.names | sort | uniq -d > refseq_finz_znf.names
./extract_gene_names.py RefSeq > ../../data/gffs/refseq_finz_znf.gff

# For Ensembl
hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/finz_seed.hmm \
    ../../data/seqs/Danio_rerio.GRCz11.pep.all.fa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/PF00096_seed.hmm \
    ../../data/seqs/Danio_rerio.GRCz11.pep.all.fa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

cat finz.names c2h2.names | sort | uniq -d > ensembl_finz_znf.names
./extract_gene_names.py Ensembl > ../../data/gffs/ensembl_finz_znf.gff
./rename_ensembl_chroms.py  ../../data/gffs/ensembl_finz_znf.gff ../../data/gffs/ensembl_finz_znf.gff

rg '>' ../../data/seqs/danio_rerio_hiqual_finz.fa | cut -c 2- | sed 's/\.t1//g' > denovo_finz_znf.names
./extract_gene_names.py denovo > ../../data/gffs/denovo_finz_znf.gff
./offset_gffs.py ../../data/gffs/denovo_finz_znf.gff ../../data/gffs/denovo_finz_znf.gff

rm *.names tmp.out

# Then calculate overlap
rg "\tmRNA\t" ../../data/gffs/refseq_finz_znf.gff \
    | gff2bed > ../../data/beds/refseq_finz_znf.transcripts.bed
rg "\tmRNA\t" ../../data/gffs/ensembl_finz_znf.gff \
    | gff2bed > ../../data/beds/ensembl_finz_znf.transcripts.bed
rg "\ttranscript\t" ../../data/gffs/denovo_finz_znf.gff \
    | gff2bed > ../../data/beds/denovo_finz_znf.transcripts.bed

bedtools intersect \
    -a ../../data/beds/ensembl_finz_znf.transcripts.bed \
    -b ../../data/beds/refseq_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.8 \
    -r \
    > ensembl_refseq.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ../../data/beds/refseq_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.8 \
    -r \
    > denovo_refseq.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ../../data/beds/ensembl_finz_znf.transcripts.bed \
    -wa \
    -wb \
    -s \
    -f 0.8 \
    -r \
    > denovo_ensembl.bed

bedtools intersect \
    -a ../../data/beds/denovo_finz_znf.transcripts.bed \
    -b ensembl_refseq.bed \
    -wa \
    -wb \
    -s \
    -f 0.8 \
    -r \
    > denovo_ensembl_refseq.bed

./count_genes.py > ../../data/finz_znf_overlap_80.txt
rm *.bed
