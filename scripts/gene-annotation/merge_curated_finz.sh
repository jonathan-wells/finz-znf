#!/usr/bin/env bash

# Purpose of this script is to combine Ensembl annotations with curated denovo
# finz-znfs and repeatmasker output. Converts to gtf file for downstream STAR

./rename_ensembl_chroms.py ../../data/gffs/Danio_rerio.GRCz11.101.gff3 subtracted.base.gff3
gffread -T subtracted.base.gff3 -o subtracted.base.gtf
gffread -T ../../data/gffs/ensembl_finz_znf.gff | awk '{ OFS="\t" }{ print $2, $3, $4, $5, $6, $7}' > subtracted.txt
rg -v -f subtracted.txt subtracted.base.gtf > subtracted.gtf
gffread -T ../../data/gffs/Danio_rerio_finz.final.gff > ../../data/gffs/Danio_rerio_finz.final.gtf 
cat subtracted.gtf ../../data/gffs/Danio_rerio_finz.final.gtf > ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf
cat ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf ../../data/gffs/danRer11.nonalt.tetranscripts.gtf > ../../data/expression/danRer11.nonalt.genes_tes.gtf

sed 's/transcript://g' ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf | sed 's/gene://g' > tmp
mv tmp ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gtf
sed 's/transcript://g' ../../data/gffs/danRer11.nonalt.tetranscripts.gtf | sed 's/gene://g' > tmp
mv tmp ../../data/gffs/danRer11.nonalt.tetranscripts.gtf
sed 's/transcript://g' ../../data/expression/danRer11.nonalt.genes_tes.gtf | sed 's/gene://g' > tmp
mv tmp ../../data/expression/danRer11.nonalt.genes_tes.gtf

rm subtracted* ../../data/gffs/Danio_rerio_finz.final.gtf

