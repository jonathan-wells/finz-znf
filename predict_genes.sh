#!/usr/bin/env bash

################################################################################
## Hefty script that predicts and collects FINZ-ZNFs.
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

for sp in ${!genomes[@]}; do
    genomefile=${genomes[$sp]}
    echo $sp
    
    # Find all chromosomes/contigs with finz domain hits.
    awk '{ 
        OFS="\t"
        } {
        if ( $11 < 1e-03 )
            if ( $9 < $10 )
                print $2, $9, $10 
            else
                print $2, $10, $9
        }' "../data/blast-out/${sp}_finz.out" |
            bedtools sort > "../data/beds/${sp}_finz.bed"

    # Extract blocks in 100kb windows around finz domains
    bedtools slop \
        -b 100000 \
        -i "../data/beds/${sp}_finz.bed" \
        -g "../data/beds/${sp}.genome" |
        bedtools merge > "../data/beds/${sp}_finz_blocks.bed"
    bedtools getfasta \
        -fi "${genomedir}/${genomefile}" \
        -bed "../data/beds/${sp}_finz_blocks.bed" > "../data/seqs/${sp}_finz_blocks.fa"

    # Predict genes
    augustus \
        --genemodel=complete \
        --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
        --proteinprofile=../data/seqs/canonical_finz_zfp_profile.prfl \
        --species=zebrafish \
         "../data/seqs/${sp}_finz_blocks.fa" > "../data/gffs/${sp}_augustus_finz.gff"

    # Extract all protein sequences
    /usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
        "../data/gffs/${sp}_augustus_finz.gff" \
        --seqfile="../data/seqs/${sp}_finz_blocks.fa"
    
    # Rename seqs by prepending species name
    for suffix in codingseq cdsexons aa; do
        mv "../data/gffs/${sp}_augustus_finz.${suffix}" "../data/seqs/${sp}_augustus_finz.${suffix}.fa"
        sed "s/>/>${sp}_/g" "../data/seqs/${sp}_augustus_finz.${suffix}.fa" > tmp.fa
        mv tmp.fa "../data/seqs/${sp}_augustus_finz.${suffix}.fa"
    done
    
done

# Extract and concatenate all FINZ-ZNFs
cat ../data/seqs/*_augustus_finz.aa.fa > tmp.fa 
hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../data/seqs/finz_domain.hmm \
    tmp.fa
awk '{ print $1 }' tmp.out | sort | uniq > finz.names
seqtk subseq tmp.fa finz.names > ../data/seqs/cypriniformes_augustus_finz.fa
rm tmp.fa tmp.out finz.names



