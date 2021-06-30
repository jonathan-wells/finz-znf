#!/usr/bin/env bash

#################################################################################
### Annotation pipeline for predicting and curating FINZ-ZNF genes.
#################################################################################

genomedir='/Users/jonwells/Genomes/Cypriniformes/'

declare -A genomes

while read line; do
    key=`awk '{ print $1 }' <<< $line`
    data=`awk '{ print $2 }' <<< $line`
    genomes[$key]="$data"
done < ../../data/species_genomes.txt

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
        -query ../../data/phmms/cypriniformes_finz_seed.consensus.fa \
        -subject "${genomedir}/${genomefile}" \
        -evalue 5e-02 \
        -outfmt 6 \
        -out "../../data/blast-out/${species}_finz_locs.out"
    
    tblastn \
        -query ../../data/phmms/cypriniformes_znf_seed.consensus.fa \
        -subject "${genomedir}/${genomefile}" \
        -evalue 5e-02 \
        -outfmt 6 \
        -out "../../data/blast-out/${species}_znf_locs.out"
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 3. Unmask regions containing finz or znf hits
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    cat "../../data/blast-out/${species}_finz_locs.out" \
        "../../data/blast-out/${species}_znf_locs.out" \
        > "../../data/blast-out/${species}_finz_znf_locs.out"
    
    # Unmask finz-znf hits with 10bp window either side.
    ./unmask.py \
        "${genomedir}/${genomefile}" \
        "../../data/blast-out/${species}_finz_znf_locs.out" \
        "${genomedir}/${species}.remasked.fa" 
    
    rm "../../data/blast-out/${species}_finz_znf_locs.out"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 4. Extract sequence from blocks >= 30kb around FINZ hits
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
        -b 15000 \
        -i "../../data/beds/${species}_finz.bed" \
        -g "../../data/beds/${species}.genome" |
        bedtools merge > "../../data/beds/${species}_finz_blocks.bed"

    # Extract sequence
    bedtools getfasta \
        -fi "${genomedir}/${species}.remasked.fa" \
        -bed "../../data/beds/${species}_finz_blocks.bed" \
        > "../../data/seqs/${species}_finz_blocks.fa"
    
    # finz-znf-unmasked genome no longer required
    rm "${genomedir}/${species}.remasked.fa"

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 5. Augustus-PPX to annotate predicted FiNZ-ZnFs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Run augustus - see config files for more detailed params
    /usr/local/Cellar/augustus/3.3.3_1/bin/augustus \
        --genemodel=complete \
        --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
        --proteinprofile=../../data/phmms/cypriniformes_finz_znf.prfl \
        --species=zebrafish \
        --minexonintronprob=0.2 \
        --softmasking=1 \
        "../../data/seqs/${species}_finz_blocks.fa" \
        > "../../data/gffs/${species}_augustus_finz.gff"
    
    # Output GFF coords correspond to finz blocks, so revert back to original
    # genome coords.
    ./offset_gffs.py \
        "../../data/gffs/${species}_augustus_finz.gff" \
        "../../data/gffs/${species}_augustus_finz.gff"

    # Extract all protein sequences
    /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
        "../../data/gffs/${species}_augustus_finz.gff" \
        --seqfile="${genomedir}/${genomefile}"
    
    # Rename seqs by prepending species name
    for suffix in "codingseq" "cdsexons" "aa"; do
        mv "../../data/gffs/${species}_augustus_finz.${suffix}" \
            "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
        sed \
            "s/>/>${species}_/g" \
            "../../data/seqs/${species}_augustus_finz.${suffix}.fa" > tmp.fa
        mv tmp.fa "../../data/seqs/${species}_augustus_finz.${suffix}.fa"
    done

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 6. HMMER to extract only predicted FINZ-ZNFs, and exclude potential
    #    pseudogenes with missing data.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    hmmsearch \
        --tblout tmp.out \
        -E 1e-02 \
        ../../data/phmms/cypriniformes_finz_seed.hmm \
        "../../data/seqs/${species}_augustus_finz.aa.fa"
    rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

    hmmsearch \
        --tblout tmp.out \
        -E 1e-02 \
        ../../data/phmms/cypriniformes_znf_seed.hmm \
        "../../data/seqs/${species}_augustus_finz.aa.fa"
    rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > znf.names

    # Generate set of labels for aa and codingseq
    cat finz.names znf.names | sort | uniq -d > finz_znf.names
   
    # Exclude genes with missing data.
    rg -B 1 'X' "../../data/seqs/${species}_augustus_finz.aa.fa" |
        rg '>' | cut -c 2- > exclude.names
    rg -v -f exclude.names finz_znf.names > tmp.names
    mv tmp.names finz_znf.names
    rm exclude.names

    # Generate labels for cdsexons
    sed "s/$/.cds/" finz_znf.names > tmp.names
    rg -f tmp.names "../../data/seqs/${species}_augustus_finz.cdsexons.fa" \
        > finz_znf.cds.names
    rm tmp.names
    
    # Filter out non-FiNZ-ZnFs from output files
    seqtk subseq \
        "../../data/seqs/${species}_augustus_finz.aa.fa" \
        finz_znf.names \
        > tmp.aa.fa
    mv tmp.aa.fa "../../data/seqs/${species}_augustus_finz.aa.fa"

    seqtk subseq \
        "../../data/seqs/${species}_augustus_finz.codingseq.fa" \
        finz_znf.names \
        > tmp.codingseq.fa
    mv tmp.codingseq.fa "../../data/seqs/${species}_augustus_finz.codingseq.fa"
    
    seqtk subseq \
        "../../data/seqs/${species}_augustus_finz.cdsexons.fa" \
        finz_znf.names \
        > tmp.cdsexons.fa
    mv tmp.cdsexons.fa "../../data/seqs/${species}_augustus_finz.cdsexons.fa"

done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. Combine species protein files and clean up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if test -f "../../data/seqs/cypriniformes_augustus_finz.fa"; then
    rm "../../data/seqs/cypriniformes_augustus_finz.fa"
fi

for species in ${!genomes[@]}; do
    cat "../../data/seqs/${species}_augustus_finz.aa.fa" \
    >> "../../data/seqs/cypriniformes_augustus_finz.fa"
done

rm tmp.out *.names

