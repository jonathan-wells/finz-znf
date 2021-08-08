#!/usr/bin/env bash

declare -A readaccs
readaccs=(
    [Danio_tinwini]=ERR3284972
)

for species in ${!readaccs[@]}; do
    echo $species
    echo ${readaccs[$species]}
    fasterq-dump ${readaccs[$species]} \
        --split-3 \
        -o "${species}_${readaccs[$species]}" \
        -e 10 \
        -p
done
        # -o "../../data/danio-reads/${species}_${readaccs[$species]}" \
