#!/usr/bin/env bash

/usr/local/Augustus/bin/augustus \
    --genemodel=complete \
    --species=zebrafish_new \
    --minexonintronprob=0.2 \
    --optCfgFile=/usr/local/Augustus/config/ppx.cfg \
    --proteinprofile=/Users/jonwells/Projects/feschottelab/finz-znf/data/phmms/drerio_finz_expressed2.prfl \
    --UTR=on \
    --softmasking=1 \
    "../../data/seqs/Danio_rerio_finz_blocks.fa" \
    > "Danio_rerio_hiqual_finz.gff"

./offset_gffs.py "Danio_rerio_hiqual_finz.gff" "Danio_rerio_hiqual_finz.gff"

/usr/local/Augustus/scripts/getAnnoFasta.pl \
    "Danio_rerio_hiqual_finz.gff" \
    --seqfile="/Users/jonwells/Genomes/Cypriniformes/GCF_000002035.6_GRCz11_genomic.nonalt.fna"

hmmsearch --tblout tmp.out \
    -E 1e-02 \
    "../../data/phmms/finz_seed.hmm" \
    "Danio_rerio_hiqual_finz.aa"
rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-02 \
    "../../data/phmms/PF00096_seed.hmm" \
    "Danio_rerio_hiqual_finz.aa"

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

cat finz.names c2h2.names | sort | uniq -d > finz_znf.names
seqtk subseq "Danio_rerio_hiqual_finz.aa" finz_znf.names > ../../data/seqs/Danio_rerio_hiqual_finz.fa
sed 's/.t1/(.t1|$)/' finz_znf.names > tmp; mv tmp finz_znf.names
rg -f finz_znf.names Danio_rerio_hiqual_finz.gff > ../../data/gffs/Danio_rerio_hiqual_finz.gff
rm c2h2.names finz.names
