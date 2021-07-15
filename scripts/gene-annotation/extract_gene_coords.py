#!/usr/bin/env python3

import offset_gffs as og
from Bio import SeqIO
import re
import sys

def get_offset_gff(species):
    filename = f'../../data/gffs/{species}_augustus_finz.gff'
    with open(filename) as infile:
        gff = [line.strip() for line in infile]
    return gff

def get_finz_genes(species):
    finz_genes = []
    for record in SeqIO.parse('../../data/seqs/cypriniformes_augustus_finz.fa', 
                              'fasta'):
        idparts = record.id.split('_')
        sp, gid = '_'.join(idparts[:2]), idparts[-1]
        if sp == species:
            finz_genes.append(gid)
    return finz_genes

def extract_gene(gff, genes):
    tpattern = '([\w\.]+)\tAUGUSTUS\ttranscript\t(\d+)\t(\d+)\t.+\tID=([\w\.]+)'
    for line in gff:
        gene = re.search(tpattern, line)
        if gene:
            chrom, start, end, gid = gene.groups()
            if gid not in genes:
                continue
            print(f'{chrom}\t{start}\t{end}\t1')

def print_circos_format(species):
    gff = get_offset_gff(species)
    finz_genes = get_finz_genes(species)
    extract_gene(gff, finz_genes)
    

if __name__ == '__main__':
    species = sys.argv[1]
    print_circos_format(species)
