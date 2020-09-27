#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import os

def load_orthos(filename):
    alignment = defaultdict(list)
    for record in SeqIO.parse(filename, 'fasta'):
        if 'X' in record.seq:
            continue
        species = '_'.join(record.id.split('_')[:2])
        alignment[species].append(record)
    for species in list(alignment.keys()):
        if len(alignment[species]) > 1:
            alignment.pop(species)
    alignment = {key: val[0] for key, val in alignment.items()}
    return alignment

def concatenate_alignments(species_list, alignments):
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

def main():
    species = ['Pimephales_promelas',
               'Cyprinus_carpio',
               'Sinocyclocheilus_anshuiensis',
               'Sinocyclocheilus_rhinocerous',
               'Sinocyclocheilus_grahami',
               'Carassius_auratus',
               'Oxygymnocypris_stewartii',
               'Anabarilius_grahami',
               'Cirrhinus_molitorella',
               'Labeo_rohita',
               'Poropuntius_huangchuchieni',
               'Hypophthalmichthys_nobilis',
               'Hypophthalmichthys_molitrix',
               'Triplophysa_siluroides',
               'Triplophysa_tibetana',
               'Culter_alburnus',
               'Megalobrama_amblycephala',
               'Onychostoma_macrolepis',
               'Labeo_catla',
               'Leuciscus_waleckii',
               'Danionella_dracula',
               'Danionella_translucida',
               'Danio_albolineatus',
               'Danio_jaintianensis',
               'Danio_choprai',
               'Danio_aesculapii',
               'Danio_kyathit',
               'Danio_tinwini',
               'Danio_rerio']
    alignments = []
    dirname = '../../data/species-phylogeny'
    for filename in os.listdir(dirname):
        if filename.endswith('trimmed.fa'):
            prot = filename.split('.')[0].split('_')[1]
            alignments.append(load_orthos(f'{dirname}/{filename}'))
    concat, meta = concatenate_alignments(species, alignments)
    with open(f'{dirname}/supermatrix_partitions.nex', 'w') as outfile:
        outfile.write(meta)
    with open(f'{dirname}/supermatrix.fa', 'w') as outfile:
        for sp in concat:
            outfile.write(f'>{sp}\n')
            outfile.write(f'{concat[sp]}\n')

if __name__ == '__main__':
    main()
