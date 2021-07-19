#!/usr/bin/env bash

/usr/local/Cellar/augustus/3.3.3_2/bin/augustus \
    --genemodel=complete \
    --species=zebrafish_new \
    --minexonintronprob=0.2 \
    --optCfgFile=/usr/local/Cellar/augustus/3.3.3_2/config/ppx.cfg \
    --proteinprofile=/Users/jonwells/Projects/feschottelab/finz-znf/data/phmms/drerio_finz_expressed2.prfl \
    --UTR=on \
    "../../data/seqs/Danio_rerio_finz_blocks.fa" \
    > "Danio_rerio_hiqual_finz.gff"

    # --softmasking=1 \
./offset_gffs.py "Danio_rerio_hiqual_finz.gff" "Danio_rerio_hiqual_finz.gff"

/usr/local/Cellar/augustus/3.3.3_2/scripts/getAnnoFasta.pl \
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
seqtk subseq "Danio_rerio_hiqual_finz.aa" finz_znf.names > ../../data/seqs/Danio_rerio_hiqual_finz.aa.fa
seqtk subseq "Danio_rerio_hiqual_finz.cdsexons" finz_znf.names > ../../data/seqs/Danio_rerio_hiqual_finz.cdsexons.fa
seqtk subseq "Danio_rerio_hiqual_finz.codingseq" finz_znf.names > ../../data/seqs/Danio_rerio_hiqual_finz.codingseq.fa
seqtk subseq "Danio_rerio_hiqual_finz.mrna" finz_znf.names >../../data/seqs/Danio_rerio_hiqual_finz.mrna.fa
sed 's/.t1/(.t1|$)/' finz_znf.names > tmp; mv tmp finz_znf.names
rg -f finz_znf.names Danio_rerio_hiqual_finz.gff > ../../data/gffs/Danio_rerio_hiqual_finz.gff

# Now replace original finz with hiqual for Danio_rerio
genomedir='/Users/jonwells/Genomes/Cypriniformes/'

declare -A genomes
# genomes[Danio_rerio]=GCF_000002035.6_GRCz11_genomic.nonalt.fna

while read line; do
    key=`awk '{ print $1 }' <<< $line`
    data=`awk '{ print $2 }' <<< $line`
    genomes[$key]="$data"
done < ../../data/species_genomes.txt

if test -f "../../data/seqs/cypriniformes_augustus_finz.fa"; then
    rm "../../data/seqs/cypriniformes_augustus_finz.fa"
fi

for species in ${!genomes[@]}; do
    cat "../../data/seqs/${species}_augustus_finz.aa.fa" \
    >> "../../data/seqs/cypriniformes_augustus_finz.fa"
done

rg '^>' ../../data/seqs/cypriniformes_augustus_finz.fa | rg -v 'Danio_rerio' | cut -c 2- > tmp.names
seqtk subseq ../../data/seqs/cypriniformes_augustus_finz.fa tmp.names > tmp.fa
sed 's/>g/>Danio_rerio_g/' ../../data/seqs/Danio_rerio_hiqual_finz.aa.fa > tmp.aa.fa
cat tmp.fa tmp.aa.fa > ../../data/seqs/cypriniformes_augustus_finz.fa

rm *.names tmp.*

