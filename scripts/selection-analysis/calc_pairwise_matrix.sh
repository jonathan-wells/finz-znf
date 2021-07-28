#!/usr/bin/env bash

seqdir="/Users/jonwells/Projects/feschottelab/finz-znf/data/seqs"
seqdump="${seqdir}/seqdump"

# Clear pre-existing seqs
rm -r $seqdump
mkdir $seqdump

# Restrict alignments to those between species with high-quality genomes.
awk '{ print $1 }' ../../data/species_genomes.txt > ../../data/hiqual_species.txt
rg -f ../../data/hiqual_species.txt "${seqdir}/cypriniformes_augustus_finz.fa" |
    cut -c 2- > hiqual.names
seqtk subseq "${seqdir}/cypriniformes_augustus_finz.fa" hiqual.names > tmp.fa
fastaexplode -f tmp.fa -d $seqdump
 
i=0
N=12
for file in $(ls $seqdump); do
    if echo $file | grep -vqe '.*\.fa' ; then
        continue
    fi
    ((i++)) 
    echo $i
    ((j=j%N)); ((j++==0)) && wait
    needle -asequence "${seqdump}/${file}" \
        -bsequence tmp.fa \
        -gapopen 12.0 \
        -gapextend 0.1 \
        -outfile "${seqdump}/needle_${i}.out" &
        
done

for file in ${seqdump}/needle*; do
    grep '# 1:' $file | awk '{ print $3 }' > tmp1.out
    grep '# 2:' $file | awk '{ print $3 }' > tmp2.out
    grep '# Identity' $file | awk '{ print $NF }' > tmp3.out
    paste -d '\t' tmp1.out tmp2.out tmp3.out >> tmp4.out
    rm tmp{1,2,3}.out
done

sed -i '.bak' 's/(//g' tmp4.out
sed -i '.bak' 's/%)//g' tmp4.out
sort -k 3 tmp4.out > ../../data/selection-analysis/pairwise_needleman.txt

################################################################################
## Cleanup
################################################################################

rm tmp*
rm hiqual.names
