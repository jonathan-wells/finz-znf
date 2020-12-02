#!/usr/bin/env bash

################################################################################
## Annotation pipeline for predicting and curating FINZ-ZNF genes.
################################################################################

chr4='/Users/jonwells/Genomes/Cypriniformes/Supplementary_dataset1_chr4.fasta'

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Extract genome contig/scaffold/chromosome lengths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bioawk -c fastx '{ 
    OFS="\t" 
} { 
    print $name, length($seq) 
}' "$chr4" > "../../data/beds/denovo_chr4.genome" 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. tblastn to get coordinates of protein domains
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tblastn \
    -query ../../data/phmms/finz_seed.consensus.fa \
    -subject "$chr4" \
    -evalue 1e-03 \
    -num_threads 8 \
    -outfmt 6 \
    -out "../../data/blast-out/denovo_chr4_finz_locs.out"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Extract sequence from blocks >= 100kb around FINZ hits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
awk '{ 
    OFS="\t"
    } {
    if ( $9 < $10 )
        print $2, $9, $10 
    else
        print $2, $10, $9
    }' "../../data/blast-out/denovo_chr4_finz_locs.out" |
        bedtools sort > "../../data/beds/denovo_chr4_finz.bed"

# Extract blocks
bedtools slop \
    -b 50000 \
    -i "../../data/beds/denovo_chr4_finz.bed" \
    -g "../../data/beds/denovo_chr4.genome" |
    bedtools merge > "../../data/beds/denovo_chr4_finz_blocks.bed"

# Extract sequence
bedtools getfasta \
    -fi "$chr4" \
    -bed "../../data/beds/denovo_chr4_finz_blocks.bed" > "../../data/seqs/denovo_chr4_finz_blocks.fa"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Augustus-ppx to annotate predicted FINZ-znfs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 augustus \
     --genemodel=complete \
     --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
     --proteinprofile=../../data/phmms/drerio_finz_expressed.prfl \
     --species=zebrafish \
      "../../data/seqs/denovo_chr4_finz_blocks.fa" > "../../data/gffs/denovo_chr4_augustus_finz.gff"

  # Extract all protein sequences
  /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
      "../../data/gffs/denovo_chr4_augustus_finz.gff" \
      --seqfile="../../data/seqs/denovo_chr4_finz_blocks.fa"

   # Rename seqs by prepending species name
   for suffix in codingseq cdsexons aa; do
       mv "../../data/gffs/denovo_chr4_augustus_finz.${suffix}" "../../data/seqs/denovo_chr4_augustus_finz.${suffix}.fa"
       sed "s/>/>denovo_chr4_/g" "../../data/seqs/denovo_chr4_augustus_finz.${suffix}.fa" > tmp.fa
       mv tmp.fa "../../data/seqs/denovo_chr4_augustus_finz.${suffix}.fa"
   done

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. HMMER to extract only predicted FINZ-ZNFs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat ../../data/seqs/denovo_chr4_augustus_finz.aa.fa > tmp.fa 
cat ../../data/seqs/denovo_chr4_augustus_finz.codingseq.fa > tmp2.fa 

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/finz_seed.hmm \
    tmp.fa
rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/PF00096_seed.hmm \
    tmp.fa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

cat finz.names c2h2.names | sort | uniq -d > finz_znf.names
seqtk subseq tmp.fa finz_znf.names > ../../data/seqs/denovo_chr4_augustus_finz.fa
seqtk subseq tmp2.fa finz_znf.names >../../data/seqs/denovo_chr4_augustus_finz.dna.fa

rm tmp*.fa tmp.out *.names

