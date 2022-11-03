#!/usr/bin/env bash

declare genomes=(
    # Anabarilius_grahami
    # Carassius_auratus
    # Cyprinus_carpio
    Danio_aesculapii
    Danio_albolineatus
    Danio_choprai
    Danio_jaintianensis
    Danio_kyathit
    Danio_rerio
    Danio_tinwini
    Danionella_dracula
    Danionella_translucida
    # Hypophthalmichthys_molitrix
    # Labeo_catla
    # Labeo_rohita
    # Leuciscus_waleckii
    # Paracanthobrama_guichenoti
    # Puntigrus_tetrazona
    # Sinocyclocheilus_anophthalmus
    # Sinocyclocheilus_anshuiensis
    # Sinocyclocheilus_rhinocerous    
)

for species in ${genomes[@]}; do
    filename="../../data/repeatmodeller-out/${species}/*_rmodeller_database-families.fa"

    if [ -f $filename ]; then
        echo $filename
        cd-hit-est \
            -i $filename \
            -c 0.8 \
            -o "../../data/repeatmodeller-out/${species}/${species}_clustered_families.fa" \
            -sf 1 \
            -sc 1 \
            -l 300 \
            -T 0 \
            -d 0 \
            -aS 0.8 \
            -aL 0.05
    fi

done
