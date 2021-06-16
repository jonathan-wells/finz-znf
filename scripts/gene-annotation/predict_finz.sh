#!/usr/bin/env bash

#################################################################################
### Annotation pipeline for predicting and curating FINZ-ZNF genes.
#################################################################################

genomedir='/Users/jonwells/Genomes/Cypriniformes/ncbi-genomes-2020-07-09'

declare -A genomes

# while read line; do
#     key=`awk '{ print $1 }' <<< $line`
#     data=`awk '{ print $2 }' <<< $line`
#     genomes[$key]="$data"
# done < ../../data/species_genomes.txt


genomes=(
    [Danio_rerio]="GCF_000002035.6_GRCz11_genomic.nonalt.remasked.fna"
)
for species in ${!genomes[@]}; do
    genomefile=${genomes[$species]}
    echo $species
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 1. Extract genome contig/scaffold/chromosome lengths
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bioawk -c fastx '{ 
        OFS="\t" 
    } { 
        print $name, length($seq) 
    }' "${genomedir}/${genomefile}" > "../../data/beds/${species}.genome" 


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 2. tblastn to get coordinates of protein domains
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tblastn \
        -query ../../data/phmms/finz_seed.consensus.fa \
        -subject "${genomedir}/${genomefile}" \
        -evalue 1e-03 \
        -outfmt 6 \
        -out "../../data/blast-out/${species}_finz_locs.out"
    
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
        }' "../../data/blast-out/${species}_finz_locs.out" |
            bedtools sort > "../../data/beds/${species}_finz.bed"

    # Extract blocks
    bedtools slop \
        -b 50000 \
        -i "../../data/beds/${species}_finz.bed" \
        -g "../../data/beds/${species}.genome" |
        bedtools merge > "../../data/beds/${species}_finz_blocks.bed"

    # Extract sequence
    bedtools getfasta \
        -fi "${genomedir}/${genomefile}" \
        -bed "../../data/beds/${species}_finz_blocks.bed" > "../../data/seqs/${species}_finz_blocks.fa"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 4. Augustus-ppx to annotate predicted FINZ-znfs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /usr/local/Cellar/augustus/3.3.3_1/bin/augustus \
        --genemodel=complete \
        --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
        --proteinprofile=../../data/phmms/drerio_finz_expressed.prfl \
        --species=zebrafish \
        --softmasking=1 \
        "../../data/seqs/${species}_finz_blocks.fa" > "../../data/gffs/${species}_augustus_finz.gff"

    # Extract all protein sequences
    /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
        "../../data/gffs/${species}_augustus_finz.gff" \
        --seqfile="../../data/seqs/${species}_finz_blocks.fa"
    
    # Rename seqs by prepending species name
    for suffix in codingseq cdsexons aa; do
        mv "../../data/gffs/${species}_augustus_finz.${suffix}" "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
        sed "s/>/>${species}_/g" "../../data/seqs/${species}_augustus_finz.${suffix}.fa" > tmp.fa
        mv tmp.fa "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
    done
done

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. HMMER to extract only predicted FINZ-ZNFs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat ../../data/seqs/*_augustus_finz.aa.fa > tmp.fa 
# cat ../../data/seqs/*_augustus_finz.codingseq.fa > tmp2.fa 

# Original is 1e-04 for both
hmmsearch --tblout tmp.out \
    -E 1e-02 \
    ../../data/phmms/finz_seed.hmm \
    tmp.fa
rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-02 \
    ../../data/phmms/PF00096_seed.hmm \
    tmp.fa
rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

# for species in ${!genomes[@]}; do
#     hmmsearch --tblout tmp.out \
#         -E 0.1 \
#         ../../data/phmms/finz_iter1.hmm \
#         "../../data/seqs/${species}_augustus_finz.aa.fa"
#     rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq >> finz.names

#     hmmsearch --tblout tmp.out \
#         -E 0.1 \
#         ../../data/phmms/PF00096_seed.hmm \
#         "../../data/seqs/${species}_augustus_finz.aa.fa"
#     rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq >> c2h2.names

# #     blastp \
# #         -query ../../data/phmms/finz_seed.consensus.fa \
# #         -subject "../../data/seqs/${species}_augustus_finz.aa.fa" \
# #         -evalue 1 \
# #         -outfmt 6 |
# #         awk '{ print $2 }' | sort | uniq >> finz.names
# #     blastp \
# #         -query ../../data/phmms/PF00096_seed.consensus.fa \
# #         -subject "../../data/seqs/${species}_augustus_finz.aa.fa" \
# #         -evalue 1 \
# #         -outfmt 6 |
# #         awk '{ print $2 }' | sort | uniq >> c2h2.names
# done

cat finz.names c2h2.names | sort | uniq -d > finz_znf.names
seqtk subseq tmp.fa finz_znf.names > ../../data/seqs/cypriniformes_augustus_finz.fa
seqtk subseq tmp2.fa finz_znf.names >../../data/seqs/cypriniformes_augustus_finz.dna.fa

rm tmp*.fa tmp.out *.names

