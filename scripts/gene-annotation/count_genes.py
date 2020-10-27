#!/usr/bin/env python3

import re

refseq = set()
ensembl = set()
denovo = 0

with open('../../data/gffs/denovo_finz_znf.gff') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            denovo += 1
    
with open('../../data/gffs/refseq_finz_znf.gff') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene = re.search(';gene=(.+?);', line[-1])
            if gene:
                refseq.add(gene.group(1))

with open('../../data/gffs/ensembl_finz_znf.gff') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene = re.search('ID=gene:(.+?);', line[-1])
            if gene:
                ensembl.add(gene.group(1))

denovo_ensembl = set()
with open('denovo_ensembl.bed') as infile:
    for line in infile:
        gene = re.search(';Parent=gene:(.+?);', line)
        if gene:
            denovo_ensembl.add(gene.group(1))

denovo_refseq = set()
with open('denovo_refseq.bed') as infile:
    for line in infile:
        gene = re.search(';gene=(.+?);', line)
        if gene:
            denovo_refseq.add(gene.group(1))

ensembl_refseq = set()
refseq_ensembl = set()
with open('ensembl_refseq.bed') as infile:
    for line in infile:
        ensembl_gene = re.search(';Parent=gene:(.+?);', line)
        refseq_gene = re.search(';Parent=gene-(.+?);', line)
        if ensembl_gene and refseq_gene:
            refseq_ensembl.add(refseq_gene.group(1))
            ensembl_refseq.add(ensembl_gene.group(1))

denovo_ensembl_refseq = set()
with open('denovo_ensembl_refseq.bed') as infile:
    for line in infile:
        ensembl_gene = re.search(';Parent=gene:(.+?);', line)
        if ensembl_gene:
            denovo_ensembl_refseq.add(ensembl_gene.group(1))

assert len(refseq_ensembl) == len(ensembl_refseq)

a = denovo - len(denovo_ensembl) - len(denovo_refseq) + len(denovo_ensembl_refseq)
b = len(ensembl) - len(denovo_ensembl) - len(ensembl_refseq) + len(denovo_ensembl_refseq)
c = len(refseq) - len(denovo_refseq) - len(ensembl_refseq) + len(denovo_ensembl_refseq)
ab = len(denovo_ensembl) - len(denovo_ensembl_refseq)
ac = len(denovo_refseq) - len(denovo_ensembl_refseq)
bc = len(ensembl_refseq) - len(denovo_ensembl_refseq)
abc = len(denovo_ensembl_refseq)

print(f'{a}\n{b}\n{ab}\n{c}\n{ac}\n{bc}\n{abc}')
