#!/usr/bin/env python

import sys
import re
from Bio import SeqIO

"""
Sub-routines for overlap_finz.sh. maps proteins to their parent genes,
enabling their extraction from gff files.
"""

def refseq_prot_to_gene():
    """Extract refseq gene ids for each protein"""
    refseqdict = {}
    with open('../../data/gffs/GCF_000002035.6_GRCz11_genomic.gff') as infile:
        for line in infile:
            protein = re.search(';protein_id=([X|N]P_\d+\.\d+)', line)
            if protein:
                gene = re.search(';gene=(.+?);', line)
                refseqdict[protein.group(1)] = gene.group(1)
    return refseqdict

def ensembl_prot_to_gene():
    """Extract ensembl gene ids for each protein"""
    ensembldict = {}
    for record in SeqIO.parse('../../data/seqs/Danio_rerio.GRCz11.pep.all.fa', 'fasta'):
        gene = re.search('gene:(ENSDARG\d+)\.', record.description).group(1)
        ensembldict[record.id] = gene
    return ensembldict

def extract_refseq_gff(proteinfile):
    refseqdict = refseq_prot_to_gene()
    with open(proteinfile) as infile:
        proteins = [line.strip() for line in infile]
    genes = set(refseqdict[i] for i in proteins)
    with open('../../data/gffs/GCF_000002035.6_GRCz11_genomic.gff') as infile:
        for line in infile:
            if line.startswith('#') and not line.startswith('##'):
                print(line.strip())
                continue
            gene = re.search(';gene=(.+?);', line)
            if gene and gene.group(1) in genes:
                print(line.strip())

def extract_ensembl_gff(proteinfile):
    ensembldict = ensembl_prot_to_gene()
    with open(proteinfile) as infile:
        proteins = [line.strip() for line in infile]
    genes = set(ensembldict[i] for i in proteins)
    ingene = False
    with open('../../data/gffs/Danio_rerio.GRCz11.101.gff3') as infile:
        print(infile.readline().strip())
        for line in infile:
            if line == '###\n':
                ingene = False
                continue
            if re.search('\t\w+\tgene\t', line):
                gene = re.search('ID=gene:(ENSDARG\d+);', line).group(1)
                if gene in genes:
                    print(line.strip())
                    ingene = True
            else:
                if ingene:
                    print(line.strip())

def extract_denovo_gff(proteinfile):
    with open(proteinfile) as infile:
        proteins = set(line.strip() for line in infile)
    with open('../../data/gffs/Danio_rerio_augustus_finz.gff') as infile:
        print(infile.readline().strip())
        for line in infile:
            if line.startswith('#'):
                continue
            
            gene = re.search('ID=(g\d+)\D', line)
            if gene:
                gene = gene.group(1)
                if f'{gene}.t1' in proteins:
                    print(line.strip())
            

if __name__ == '__main__':
    if sys.argv[1] == 'RefSeq':
        extract_refseq_gff('refseq_finz_znf.names')
    elif sys.argv[1] == 'Ensembl':
        extract_ensembl_gff('ensembl_finz_znf.names')
    elif sys.argv[1] == 'denovo':
        extract_denovo_gff('denovo_finz_znf.names')
