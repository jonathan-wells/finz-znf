#!/usr/bin/env bash

seqdir="/Users/jonwells/Projects/feschottelab/finz-znf/data/seqs"
seqdump="${seqdir}/seqdump_denovo_chr4"

# fastaexplode -f "${seqdir}/denovo_chr4_augustus_finz.fa" -d $seqdump
 
# i=0
# for file in `ls $seqdump`; do
#     if echo $file | grep -vqe '.*\.fa' ; then
#         continue
#     fi
#     echo $i
#     needle -asequence "${seqdump}/${file}" \
#         -bsequence "${seqdir}/denovo_chr4_augustus_finz.fa" \
#         -gapopen 12.0 \
#         -gapextend 0.1 \
#         -outfile "${seqdump}/denovo_chr4_needle${i}.out"
        
#     ((i=i+1))
# done

for file in ${seqdump}/denovo_chr4_needle*; do
    # if echo $file | rg -v "denovo_chr4_needle.*\.out" ; then
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
sort -k 3 tmp4.out > ../../data/denovo_chr4_pairwise_needleman.txt

################################################################################
## Cleanup
################################################################################

rm tmp*

