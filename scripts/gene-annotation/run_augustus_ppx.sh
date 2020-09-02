#!/usr/bin/env bash

################################################################################
## Annotation pipeline for predicting and curating FINZ-ZNF genes.
################################################################################

genomedir='/Users/jonwells/Genomes/Cypriniformes/ncbi-genomes-2020-07-09'

declare -A genomes

while read line; do
    key=`awk '{ print $1 }' <<< $line`
    data=`awk '{ print $2 }' <<< $line`
    genomes[$key]="$data"
done < ../../data/species_genomes.txt

#for species in ${!genomes[@]}; do
#    genomefile=${genomes[$species]}
#    echo $species
    
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    # 1. Extract genome contig/scaffold/chromosome lengths
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    bioawk -c fastx '{ 
#        OFS="\t" 
#    } { 
#        print $name, length($seq) 
#    }' "${genomedir}/${genomefile}" > "../../data/beds/${species}.genome" 


#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    # 2. tblastn to get coordinates of protein domains
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    tblastn \
#        -query ../../data/phmms/finz_seed.consensus.fa \
#        -subject "${genomedir}/${genomefile}" \
#        -evalue 1e-03 \
#        -num_threads 8 \
#        -outfmt 6 \
#        -out "../../data/blast-out/${species}_finz_locs.out"

#    for prot in FICD GLYT ILF2 RAG1 RAG2 TBR1; do
#        tblastn \
#            -query "../../data/phmms/${prot}.consensus.fa" \
#            -subject "${genomedir}/${genomefile}" \
#            -evalue 1e-05 \
#            -num_threads 8 \
#            -outfmt 6 \
#            -out "../../data/blast-out/${species}_${prot}_locs.out"
#    done

#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    # 3. Extract sequence from blocks >= 100kb around FINZ hits
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    awk '{ 
#        OFS="\t"
#        } {
#        if ( $9 < $10 )
#            print $2, $9, $10 
#        else
#            print $2, $10, $9
#        }' "../../data/blast-out/${species}_finz_locs.out" |
#            bedtools sort > "../../data/beds/${species}_finz.bed"

#    # Extract blocks
#    bedtools slop \
#        -b 50000 \
#        -i "../../data/beds/${species}_finz.bed" \
#        -g "../../data/beds/${species}.genome" |
#        bedtools merge > "../../data/beds/${species}_finz_blocks.bed"

#    # Extract sequence
#    bedtools getfasta \
#        -fi "${genomedir}/${genomefile}" \
#        -bed "../../data/beds/${species}_finz_blocks.bed" > "../../data/seqs/${species}_finz_blocks.fa"
    
#    for prot in FICD GLYT ILF2 RAG1 RAG2 TBR1; do
#        awk '{ 
#            OFS="\t"
#            } {
#            if ( $9 < $10 )
#                print $2, $9, $10 
#            else
#                print $2, $10, $9
#            }' "../../data/blast-out/${species}_${prot}_locs.out" |
#                bedtools sort > "../../data/beds/${species}_${prot}.bed"

#        # Extract blocks
#        bedtools slop \
#            -b 50000 \
#            -i "../../data/beds/${species}_${prot}.bed" \
#            -g "../../data/beds/${species}.genome" |
#            bedtools merge > "../../data/beds/${species}_${prot}_blocks.bed"

#        # Extract sequence
#        bedtools getfasta \
#            -fi "${genomedir}/${genomefile}" \
#            -bed "../../data/beds/${species}_${prot}_blocks.bed" \
#            > "../../data/seqs/${species}_${prot}_blocks.fa"
#    done

#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    # 4. Augustus-ppx to annotate predicted FINZ-znfs
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     augustus \
#         --genemodel=complete \
#         --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
#         --proteinprofile=../../data/phmms/drerio_finz_expressed.prfl \
#         --species=zebrafish \
#          "../../data/seqs/${species}_finz_blocks.fa" > "../../data/gffs/${species}_augustus_finz.gff"

#      Extract all protein sequences
#      /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
#          "../../data/gffs/${species}_augustus_finz.gff" \
#          --seqfile="../../data/seqs/${species}_finz_blocks.fa"
    
#       # Rename seqs by prepending species name
#       for suffix in codingseq cdsexons aa; do
#           mv "../../data/gffs/${species}_augustus_finz.${suffix}" "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
#           sed "s/>/>${species}_/g" "../../data/seqs/${species}_augustus_finz.${suffix}.fa" > tmp.fa
#           mv tmp.fa "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
#       done
    
#    for prot in FICD GLYT ILF2 RAG1 RAG2 TBR1; do
#        augustus \
#            --genemodel=complete \
#            --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
#            --proteinprofile="../../data/phmms/${prot}.prfl" \
#            --species=zebrafish \
#             "../../data/seqs/${species}_${prot}_blocks.fa" \
#             > "../../data/gffs/${species}_augustus_${prot}.gff"

#         # Extract all protein sequences
#         /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
#             "../../data/gffs/${species}_augustus_${prot}.gff" \
#             --seqfile="../../data/seqs/${species}_${prot}_blocks.fa"
        
#          # Rename seqs by prepending species name
#          for suffix in codingseq cdsexons aa; do
#              mv "../../data/gffs/${species}_augustus_${prot}.${suffix}" "../../data/seqs/${species}_augustus_${prot}.${suffix}.fa"
#              sed "s/>/>${species}_/g" "../../data/seqs/${species}_augustus_${prot}.${suffix}.fa" > tmp.fa
#              mv tmp.fa "../../data/seqs/${species}_augustus_${prot}.${suffix}.fa"
#          done

#    done

#done

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. HMMER to extract only predicted FINZ-ZNFs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#cat ../../data/seqs/*_augustus_finz.aa.fa > tmp.fa 
#cat ../../data/seqs/*_augustus_finz.codingseq.fa > tmp2.fa 

#hmmsearch --tblout tmp.out \
#    -E 1e-04 \
#    ../../data/phmms/finz_seed.hmm \
#    tmp.fa
#rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

#hmmsearch --tblout tmp.out \
#    -E 1e-04 \
#    ../../data/phmms/PF00096_seed.hmm \
#    tmp.fa

#rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

#cat finz.names c2h2.names | sort | uniq -d > finz_znf.names
#seqtk subseq tmp.fa finz_znf.names > ../../data/seqs/cypriniformes_augustus_finz.fa
#seqtk subseq tmp2.fa finz_znf.names >../../data/seqs/cypriniformes_augustus_finz.codingseq.fa

# Use blastp to extract only homologues of marker genes and discard flanking
# genes captured by augustus.
for prot in FICD GLYT ILF2 RAG1 RAG2 TBR1; do
    cat ../../data/seqs/*_augustus_$prot.aa.fa \
        > "../../data/species-phylogeny/multispecies_${prot}.fa"
    blastp \
        -query "../../data/phmms/${prot}.consensus.fa" \
        -subject "../../data/species-phylogeny/multispecies_${prot}.fa" \
        -outfmt 6 \
        -evalue 1e-10 | awk '{ print $2 }' > "${prot}.names"
    seqtk subseq \
        "../../data/species-phylogeny/multispecies_${prot}.fa" \
        "${prot}.names" > tmp.fa
    mv tmp.fa "../../data/species-phylogeny/multispecies_${prot}.fa"
done

rm tmp.out *.names

