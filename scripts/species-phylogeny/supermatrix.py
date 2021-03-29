#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os

def load_orthos(filename, species_list, mintaxa):
    alignment = defaultdict(list)
    for record in SeqIO.parse(filename, 'fasta'):
        species = '_'.join(record.id.split('_')[:2])
        if species not in species_list:
            continue
        alignment[species].append(record)
    for species in list(alignment.keys()):
        if len(alignment[species]) > 1:
            alignment.pop(species)
    alignment = {key: val[0] for key, val in alignment.items()}
    if len(alignment) < mintaxa:
        return None
    return alignment

def concatenate_alignments(species_list, alignments):
    """Produces a concatenated alignment for iqtree."""
    concatenated = defaultdict(str)
    metadata = ['#nexus', 'begin sets;']
    currpos = 0
    i = 1
    for prot in alignments:
        length = len(list(prot.values())[0].seq)
        metadata.append(f'\tcharset part{i} = {currpos+1}-{currpos+length};')
        i += 1
        currpos += length
        for species in species_list:
            record = prot.get(species, False)
            if record:
                seq = record.seq
            else:
                seq = '-'*length
            concatenated[species] += seq
    metadata.append('end;')
    metadata = '\n'.join(metadata)
    return concatenated, metadata

def clean_alignments(species_list, alignments, min_species):
    species_set = set(species_list)
    cleaned_alignments = defaultdict(list)
    for ali in alignments:
        ali_species = set(ali.keys())
        length = len(list(ali.values())[0].seq)
        busco = list(ali.values())[0].name.split('_')[-1]
        nspecies = len(species_set.intersection(ali_species))
        if nspecies < min_species:
            continue
        
        for sp in species_list:
            record = ali.get(sp, False)
            if record:
                record.id = sp
                record.name = sp
                record.description = busco
            else:
                record = SeqRecord(
                               Seq('-'*length),
                               id=sp,
                               name=sp,
                               description=busco
                               )
            cleaned_alignments[busco].append(record)
    return cleaned_alignments

def main():
    with open('../../data/species_genomes.txt') as infile:
        species = [line.split()[0] for line in infile]
    
    # Exclude these species because of poor genome quality
    exclude = ['Pimephales_promelas',
               'Cirrhinus_molitorella',
               'Labeo_gonius',
               'Poropuntius_huangchuchieni',
               'Hypophthalmichthys_nobilis',
               'Paedocypris_carbunculus',
               'Paedocypris_micromegethes']
    for sp in exclude:
        species.remove(sp)
    
    alignments = []
    dirname = '../../data/species-phylogeny'
    for filename in os.listdir(f'{dirname}/aligned-busco'):
        if filename.endswith('trimmed.fa'):
            alignment = load_orthos(f'{dirname}/aligned-busco/{filename}',
                                    species, 
                                    25)
            if alignment != None:
                alignments.append(alignment)
    
    # Write iqtree data
    concat, meta = concatenate_alignments(species, alignments)
    with open(f'{dirname}/iqtree-data/busco_supermatrix_partitions.nex', 'w') as outfile:
        outfile.write(meta)
    with open(f'{dirname}/iqtree-data/busco_supermatrix.fa', 'w') as outfile:
        for sp in concat:
            outfile.write(f'>{sp}\n')
            outfile.write(f'{concat[sp]}\n')
    
    # Write RevBayes data
    min_species = 25
    cleaned_alignments = clean_alignments(species, alignments, min_species)
    for busco, alignment in cleaned_alignments.items():
        SeqIO.write(alignment, 
                    f'{dirname}/revbayes-data/busco-input/{busco}.fa', 
                    'fasta')
    

if __name__ == '__main__':
    main()
