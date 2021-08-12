#!/usr/bin/env bash

################################################################################
## Genome files
################################################################################

genomedir="/workdir/Genomes/Cypriniformes"
declare -A genomes=(
    ["Puntigrus_tetrazona"]="GCA_004120215.1_ASM412021v1_genomic.fna"
    ["Labeo_rohita"]="GCA_018831695.1_ASM1883169v1_genomic.fna"
    ["Paracanthobrama_guichenoti"]="GCA_018749465.1_ASM1874946v1_genomic.fna"
    # "fDanAes4.curated_primary.20190724.fa"
    # "fDanDra1.1.fa"
    # "fDanKya3.curated_primary.20190802.fa"
    # "fDanTin1.curated_primary.20190226.fa"
    # "fDanJai1.curated_primary.20190227.fa"
    # "fDanTra1.curated_assembly.fa"
    # "fDanCho1.curated_primary.20190226.fa"
    # "fDanAlb1.curated_primary.20190208.fa"
    # "fDreABz2.toplevel_curated.20190410.fa"
    # "fDreCBz1.curated_primary.toplevel.20190604.fa"
    # "fDreNAz1.curated_primary.toplevel.20190604.fa"
    )

for species in ${!genomes[@]}; do
    
    ############################################################################
    ## Create output directories
    ############################################################################
    genomefile=${genomes[$species]}
    
    if [ ! -d "repeatmodeller-out" ]; then
        mkdir "repeatmodeller-out"
    fi
    
    if [ ! -d "repeatmodeller-out/${species}" ]; then 
        mkdir "repeatmodeller-out/${species}"
    fi
    
    ############################################################################
    ## Run RepeatModeller
    ############################################################################

    echo "RepeatModeller - processing ${species}"
    
    /programs/RepeatModeler-2.0.1/BuildDatabase \
    -name "repeatmodeller-out/${species}/${species}_rmodeller_database" \
    "${genomedir}/${genomefile}"

    /programs/RepeatModeler-2.0.1/RepeatModeler \
    -database "repeatmodeller-out/${species}/${species}_rmodeller_database" \
    -LTRStruct \
    -pa 16
    
    echo "RepeatModeller - done ${species}"

done

