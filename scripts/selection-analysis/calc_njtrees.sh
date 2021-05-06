#!/usr/bin/env bash

# distmats=(
#     D2k4
#     D2nk4
#     D2sk4
#     D2stk4
#     D2k8
#     D2nk8
#     D2sk8
#     D2stk8
# )

# for i in ${distmats[@]}; do
#     fneighbor \
#         -datafile "../../data/selection-analysis/cypriniformes_augustus_finz_${i}.mat" \
#         -outfile "../../data/selection-analysis/cypriniformes_augustus_finz_${i}" \
#         -matrixtype l \
#         -treetype n &
# done

fneighbor \
    -datafile "../../data/selection-analysis/pairwise_needleman.mat" \
    -outfile "../../data/selection-analysis/pairwise_needleman.nwk" \
    -matrixtype l \
    -treetype n &
