#!/usr/bin/env bash

seqdir="/Users/jonwells/Projects/feschottelab/finz-znf/data/seqs"
seqdump="${seqdir}/seqdump"

# fastaexplode -f "${seqdir}/cypriniformes_augustus_finz.fa" -d $seqdump
 
# i=0
# for file in `ls $seqdump`; do
#     if echo $file | grep -vqe '.*\.fa' ; then
#         continue
#     fi
#     echo $i
#     water -asequence "${seqdump}/${file}" \
#         -bsequence "${seqdir}/cypriniformes_augustus_finz.fa" \
#         -gapopen 12.0 \
#         -gapextend 0.0 \
#         -outfile "${seqdump}/water${i}.out"
        
#     ((i=i+1))
# done

for file in $seqdump/water*; do
    # if echo $file | grep -vqe 'water.*\.out' ; then
    #     continue
    # fi
    grep '# 1:' $file | awk '{ print $3 }' > tmp1.out
    grep '# 2:' $file | awk '{ print $3 }' > tmp2.out
    grep '# Identity' $file | awk '{ print $NF }' > tmp3.out
    paste -d '\t' tmp1.out tmp2.out tmp3.out >> tmp4.out
    rm tmp{1,2,3}.out
done

sed -i '.bak' 's/(//g' tmp4.out
sed -i '.bak' 's/%)//g' tmp4.out
sort -k 3 tmp4.out > ../../data/pairwise_waterman.txt

################################################################################
## Cleanup
################################################################################

rm tmp*

