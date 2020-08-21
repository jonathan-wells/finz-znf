#!/usr/bin/env bash

################################################################################
## Annotation pipeline for predicting and curating FINZ-ZNF genes.
################################################################################

genomedir='/Users/jonwells/Genomes/Cypriniformes/ncbi-genomes-2020-07-09'
declare -A genomes=(
    [Pimephales_promelas]=GCA_000700825.1_FHM_SOAPdenovo_genomic.fna
    [Cyprinus_carpio]=GCA_000951615.2_common_carp_genome_genomic.fna
    [Sinocyclocheilus_anshuiensis]=GCA_001515605.1_SAMN03320099.WGS_v1.1_genomic.fna
    [Sinocyclocheilus_rhinocerous]=GCA_001515625.1_SAMN03320098_v1.1_genomic.fna
    [Sinocyclocheilus_grahami]=GCA_001515645.1_SAMN03320097.WGS_v1.1_genomic.fna
    [Carassius_auratus]=GCA_003368295.1_ASM336829v1_genomic.fna
    [Oxygymnocypris_stewartii]=GCA_003573665.1_Novo_Ost_1.0_genomic.fna
    [Anabarilius_grahami]=GCA_003731715.1_BGI_Agra_1.0_genomic.fna
    [Cirrhinus_molitorella]=GCA_004028445.1_ASM402844v1_genomic.fna
    [Labeo_rohita]=GCA_004120215.1_ASM412021v1_genomic.fna
    [Poropuntius_huangchuchieni]=GCA_004124795.1_ASM412479v1_genomic.fna
    [Hypophthalmichthys_nobilis]=GCA_004193235.1_HypNob1.0_genomic.fna
    [Hypophthalmichthys_molitrix]=GCA_004764525.1_HypMol1.0_genomic.fna
    [Triplophysa_siluroides]=GCA_006030095.1_ASM603009v1_genomic.fna
    [Triplophysa_tibetana]=GCA_008369825.1_ASM836982v1_genomic.fna
    [Culter_alburnus]=GCA_009869775.1_ASM986977v1_genomic.fna
    [Megalobrama_amblycephala]=GCA_009869865.1_ASM986986v1_genomic.fna
    [Onychostoma_macrolepis]=GCA_012432095.1_ASM1243209v1_genomic.fna
    [Labeo_catla]=GCA_012976165.1_CIFA_Catla_1.0_genomic.fna
    [Leuciscus_waleckii]=GCA_900092035.1_Amur_ide_genome_genomic.fna
    [Danionella_dracula]=GCA_900490495.1_fDanDra1.1_genomic.fna
    [Danionella_translucida]=GCA_903798025.1_fDanTra1.1_genomic.fna
    [Danio_albolineatus]=GCA_903798035.1_fDanAlb1.1_genomic.fna
    [Danio_jaintianensis]=GCA_903798115.1_fDanJai1.1_genomic.fna
    [Danio_choprai]=GCA_903798125.1_fDanCho1.1_genomic.fna
    [Danio_aesculapii]=GCA_903798145.1_fDanAes4.1_genomic.fna
    [Danio_kyathit]=GCA_903798195.1_fDanKya3.1_genomic.fna
    [Danio_tinwini]=GCA_903798205.1_fDanTin1.1_genomic.fna
    [Danio_rerio]=GCF_000002035.6_GRCz11_genomic.fna
)


for species in ${!genomes[@]}; do
    genomefile=${genomes[$species]}
    #echo $species
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 1. Extract genome contig/scaffold/chromosome lengths
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #bioawk -c fastx '{ 
    #    OFS="\t" 
    #} { 
    #    print $name, length($seq) 
    #}' "${genomedir}/${genomefile}" > "../../data/beds/${species}.genome" 


    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 2. tblastn to get coordinates of FINZ domains
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #tblastn \
    #    -query ../../data/phmms/finz_seed.consensus.fa \
    #    -subject "${genomedir}/${genomefile}" \
    #    -outfmt 6 \
    #    -out "../../data/blast-out/${species}_finz_locs.out"

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 3. Extract sequence from blocks >= 100kb around FINZ hits
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #awk '{ 
    #    OFS="\t"
    #    } {
    #    if ( $11 < 1e-03 )
    #        if ( $9 < $10 )
    #            print $2, $9, $10 
    #        else
    #            print $2, $10, $9
    #    }' "../../data/blast-out/${species}_finz_locs.out" |
    #        bedtools sort > "../../data/beds/${species}_finz.bed"

    ## Extract blocks
    #bedtools slop \
    #    -b 50000 \
    #    -i "../../data/beds/${species}_finz.bed" \
    #    -g "../../data/beds/${species}.genome" |
    #    bedtools merge > "../../data/beds/${species}_finz_blocks.bed"

    ## Extract sequence
    #bedtools getfasta \
    #    -fi "${genomedir}/${genomefile}" \
    #    -bed "../../data/beds/${species}_finz_blocks.bed" > "../../data/seqs/${species}_finz_blocks.fa"

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## 4. Augustus-ppx to annotate predicted FINZ-znfs
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #augustus \
    #    --genemodel=complete \
    #    --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
    #    --proteinprofile=../../data/phmms/drerio_finz_expressed.prfl \
    #    --species=zebrafish \
    #     "../../data/seqs/${species}_finz_blocks.fa" > "../../data/gffs/${species}_augustus_finz.gff"

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. HMMER to extract only predicted FINZ-ZNFs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat ../../data/seqs/*_augustus_finz.aa.fa > tmp.fa 
cat ../../data/seqs/*_augustus_finz.codingseq.fa > tmp2.fa 

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
seqtk subseq tmp.fa finz_znf.names > ../../data/seqs/cypriniformes_augustus_finz.fa
seqtk subseq tmp2.fa finz_znf.names >../../data/seqs/cypriniformes_augustus_finz.codingseq.fa

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Remove Danio rerio alternate haplotype-derived sequences
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rg -f ../../data/misc/Danio_rerio_alt_contigs.txt ../../data/gffs/Danio_rerio_augustus_finz.gff |
    rg transcript | 
    awk '{ print $NF }' | 
    sed 's/^ID=/Danio_rerio_/g' |
    sed 's/;Parent=.*//g' > tmp.out
rg '>' ../../data/seqs/cypriniformes_augustus_finz.fa |
    cut -c 2- |
    rg -v -f tmp.out > finz.names
seqtk subseq ../../data/seqs/cypriniformes_augustus_finz.fa finz.names > tmp.fa
mv tmp.fa ../../data/seqs/cypriniformes_augustus_finz.fa
seqtk subseq ../../data/seqs/cypriniformes_augustus_finz.codingseq.fa finz.names > tmp.fa
mv tmp.fa ../../data/seqs/cypriniformes_augustus_finz.codingseq.fa

rm tmp.out *.names

