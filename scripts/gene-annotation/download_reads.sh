#!/usr/bin/env bash

declare -A readaccs
readaccs=(
    [Danio_albolineatus]=ERR3284960
    [Danio_jaintianensis]=ERR3284968
    [Danio_choprai]=ERR3284956
    [Danio_tinwini]=ERR3284972
)

for species in ${!readaccs[@]}; do
    echo $species
    echo ${readaccs[$species]}
    fasterq-dump ${readaccs[$species]} \
        --split-3 \
        -o "../../data/danio-reads/${species}_${readaccs[$species]}" \
        -e 10 \
        -p
done
