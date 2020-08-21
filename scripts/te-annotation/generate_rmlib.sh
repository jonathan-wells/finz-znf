#!/usr/bin/env bash

declare -A species=(
    [Danionella_dracula]=DanDra1
    [Danionella_translucida]=DanTra1
    [Danio_albolineatus]=DanAlb1
    [Danio_jaintianensis]=DanJai1
    [Danio_choprai]=DanCho1
    [Danio_aesculapii]=DanAes4
    [Danio_kyathit]=DanKya3
    [Danio_tinwini]=DanTin1
)

# Cheeky reuse of old repbase libraries.
/programs/RepeatMasker/util/buildRMLibFromEMBL.pl \
    /programs/RepeatMasker_4-0-8/Libraries/RepeatMaskerLib.embl \
    > ../../data/rmlibs/rmlib_plus_danioninae.fa

# Add on de-novo RepeatModeller Danioninae libs
for spec in ${!species[@]}; do
    rmlib="../../data/rmlibs/${species[$spec]}_rmodeller_database-families.fa"
    ./rename_rmodeller_families.py $spec $rmlib \
        >> ../../data/rmlibs/rmlib_plus_danioninae.fa
done

