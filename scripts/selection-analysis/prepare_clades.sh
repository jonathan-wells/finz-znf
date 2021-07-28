#!/usr/bin/env bash

./extract_clades.py
nclades=$(ls ../../data/selection-analysis/pairwise_needleman*.taxa | wc -l )
for i in $(seq 1 $nclades); do
    cdir="../../data/selection-analysis/clade${i}"
    
    # Create directories and clean gene labels
    mkdir $cdir
    mv "../../data/selection-analysis/pairwise_needleman_clade${i}.taxa" "$cdir"
    sed 's/_\[.*$//' "${cdir}/pairwise_needleman_clade${i}.taxa" | 
        sed 's/Danio_rerio_//' > "${cdir}/tmp.names"
    
    # Extract seqs
    seqtk subseq \
        "../../data/seqs/Danio_rerio_hiqual_finz.codingseq.fa" \
        "${cdir}/tmp.names" > "${cdir}/pairwise_needleman_clade${i}.fa"
    
    # Align, trim, and create tree.
    prank \
        -d="${cdir}/pairwise_needleman_clade${i}.fa" \
        -o="${cdir}/pairwise_needleman_clade${i}.prank" \
        -codon
    mv "${cdir}/pairwise_needleman_clade${i}.prank.best.fas" "${cdir}/pairwise_needleman_clade${i}.prank.fa"
    trimal \
        -in "${cdir}/pairwise_needleman_clade${i}.prank.fa" \
        -out "${cdir}/pairwise_needleman_clade${i}.final.fa" \
        -nogaps
    iqtree2 \
        -s "${cdir}/pairwise_needleman_clade${i}.final.fa" \
        -st CODON

    # Prepare clade-specific template files. These use values of kappa and omega
    # from iqtree output
    kappa=$(rg "\(kappa\): ([0-9\.]+)" -or '$1' "${cdir}/pairwise_needleman_clade${i}.final.fa.log" | tail -n 1)
    omega=$(rg "\(omega\): ([0-9\.]+)" -or '$1' "${cdir}/pairwise_needleman_clade${i}.final.fa.log" | tail -n 1)
    sed "s/seqfile =/seqfile = pairwise_needleman_clade${i}.final.fa/" codeml.ctl > tmp1.ctl
    sed "s/treefile =/treefile = pairwise_needleman_clade${i}.final.fa.treefile/" tmp1.ctl > tmp2.ctl
    sed "s/outfile =/outfile = pairwise_needleman_clade${i}.final.out/" tmp2.ctl > tmp3.ctl
    sed "s/ kappa =/ kappa = $kappa/" tmp3.ctl > tmp4.ctl
    sed "s/ omega =/ omega = $omega/" tmp4.ctl > "${cdir}/codeml.ctl"
    
    # Cleanup
    rm tmp*.ctl
    rm ${cdir}/tmp.names

done
