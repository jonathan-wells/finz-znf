#!/usr/bin/env python3

import os
import re
from Bio import SeqIO

def get_species_names():
    species_dict = {}
    with open('../../data/species-phylogeny/Results_Sep08/WorkingDirectory/SequenceIDs.txt') as infile:
        for line in infile:
            seqid = re.match('\d_\d+:\s(.+?)\s', line).group(1)
            species = re.search('\[(\w+\s\w+)\]$', line).group(1)
            species = species.replace(' ', '_')
            species_dict[seqid] = species
    return species_dict

def get_gene_name():
    with open('../../data/species-phylogeny/orthogroup-consensi/og_annotations.txt') as infile:
        gene_names = {line.split()[1].upper(): line.split()[0] for line in infile}
    return gene_names

def load_single_copy_orthologues():
    with open('../../data/species-phylogeny/Results_Sep08/Single_Copy_Orthologue_Sequences') as infile:
        pass

def main():
    gd = get_gene_name()
    sd = get_species_names()
    for seq in ['TPO', 'FSHR', 'WFIKKN1']:
        og = gd[seq]
        with open(f'../../data/species-phylogeny/of_multispecies_{seq}.fa', 'w') as outfile:
            for record in SeqIO.parse(f'../../data/species-phylogeny/Results_Sep08/Single_Copy_Orthologue_Sequences/{og}.fa', 'fasta'):
                record.id = sd[record.id]
                outfile.write(f'>{record.id}\n{str(record.seq)}\n')

if __name__ == '__main__':
    main()
