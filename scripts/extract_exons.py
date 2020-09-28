#!/usr/bin/env python3

from Bio import SeqIO
# from Bio.Seq import Seq
from collections import defaultdict, deque
import re

def get_finz():
    finz = defaultdict(list)
    for record in SeqIO.parse('../data/seqs/cypriniformes_augustus_finz.fa', 'fasta'):
        seqid = record.id.split('_')
        species, gene = '_'.join(seqid[:2]), seqid[2]
        finz[species].append(gene)
    return finz

def load_cds(species, genes):
    cds = defaultdict(deque)
    with open(f'../data/gffs/{species}_augustus_finz.gff') as infile:
        for line in infile:
            if re.match('^#', line):
                continue
            line = line.strip().split('\t')
            if line[2] != 'CDS':
                continue
            parent = re.search('ID=.+\.cds;Parent=(.+)', line[-1]).group(1)
            if parent not in genes:
                continue
            strand = line[6]
            if strand == '+':
                cds[(parent, line[0])].append((int(line[3]), int(line[4])))
            else:
                cds[(parent, line[0])].appendleft((int(line[4]), int(line[3])))
    return cds

def extract_exons(species, cds):
    exons = defaultdict(list)
    blocks = SeqIO.to_dict(SeqIO.parse(f'../data/seqs/{species}_finz_blocks.fa', 'fasta'))
    for gene, exon_coords in cds.items():
        gene, block = gene
        print(f'>{species}_{gene}')
        for i, j in exon_coords:
            # print(i, j)
            if j > i:
                print(blocks[block].seq[i-1:j])
                # print(blocks[block].seq[i-1:j].translate())
                break
            else:
                # print(blocks[block].seq[j-1:i].reverse_complement().translate())
                print(blocks[block].seq[j-1:i].reverse_complement())
                break
                # print(blocks[block].seq[j-1:i].translate())

        

def main():
    with open('../data/species_genomes.txt') as infile:
        species_list = [line.split()[0] for line in infile]
    for species in species_list:
        finz = get_finz()
        cds = load_cds(species, finz[species])
        extract_exons(species, cds)

if __name__ == '__main__':
    main()
