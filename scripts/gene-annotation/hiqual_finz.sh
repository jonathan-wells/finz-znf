#!/usr/bin/env bash

blat \
    -minIdentity=99 \
    ../../data/seqs/Danio_rerio_finz_blocks.fa \
    ../../data/expression/expressed.cds.fa \
    expressed.cds.psl
pslCDnaFilter -maxAligns=1 expressed.cds.psl expressed.cds.f.psl
/usr/local/Cellar/augustus/3.3.3_1/scripts/blat2hints.pl --in=expressed.cds.f.psl --out=hints.E.gff

augustus \
    --genemodel=complete \
    --optCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/ppx.cfg \
    --extrinsicCfgFile=/usr/local/Cellar/augustus/3.3.3_1/config/extrinsic/extrinsic.MPE.cfg \
    --proteinprofile=../../data/phmms/drerio_finz_expressed.prfl \
    --hintsfile=../../data/expression/hints.E.gff \
    --species=zebrafish \
    --codingseq=on \
    --protein=on \
    --outfile=danio_rerio_hiqual_augustus_finz.gff \
    ../../data/seqs/Danio_rerio_finz_blocks.fa

/usr/local/Cellar/augustus/3.3.3_1/scripts/getAnnoFasta.pl \
    "danio_rerio_hiqual_augustus_finz.gff" \
    --seqfile="../../data/seqs/Danio_rerio_finz_blocks.fa"

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/finz_seed.hmm \
    danio_rerio_hiqual_augustus_finz.aa
rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > finz.names

hmmsearch --tblout tmp.out \
    -E 1e-04 \
    ../../data/phmms/PF00096_seed.hmm \
    danio_rerio_hiqual_augustus_finz.aa

rg -v '^#' tmp.out | awk '{ print $1 }' | sort | uniq > c2h2.names

cat finz.names c2h2.names | sort | uniq -d > finz_znf.names
seqtk subseq danio_rerio_hiqual_augustus_finz.aa finz_znf.names > danio_rerio_finz.fa
