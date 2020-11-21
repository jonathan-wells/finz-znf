#!/usr/bin/env bash

# Purpose of this script is to combine Ensembl gff files with curated denovo
# finz-znf genes.

./rename_ensembl_chroms.py ../../data/gffs/Danio_rerio.GRCz11.101.gff3 subtracted.base.gff3
awk '{ OFS="\t" }{ print $2, $3, $4, $5, $6, $7}' ../../data/gffs/ensembl_finz_znf.gff > subtracted.txt
rg -v -f subtracted.txt subtracted.base.gff3 > subtracted.gff3
cat subtracted.gff3 ../../data/gffs/Danio_rerio_finz.final.gff > subtracted.merged.gff3
bedtools sort -header -i subtracted.merged.gff3 > ../../data/gffs/Danio_rerio.GRCz11.101.curated_finz.gff3
rm subtracted*

