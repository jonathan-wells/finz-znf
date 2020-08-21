#!/usr/bin/env python3

import sys
from Bio import SeqIO

# Generate unique initial suffix for RepeatModeller consensi
species = sys.argv[1].replace('_', ' ')
initials = species.split()
initials = f'{initials[0][0].upper()}{initials[1][0].upper()}{initials[1][1:3]}'

# Rename them seqs
for record in SeqIO.parse(sys.argv[2], 'fasta'):
    id_parts = record.id.split('#')
    record.id = f'{id_parts[0]}_{initials}#{id_parts[1]}' 
    record.description = f'@{species}  [S:] RepbaseID: XX'
    SeqIO.write(record, sys.stdout, 'fasta')

