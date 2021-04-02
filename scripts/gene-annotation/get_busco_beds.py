#!/usr/bin/env python

import os
from Bio import SeqIO
import re

scbuscos = {}
for species in ['Danio_albolineatus',
                'Danio_choprai',
                'Danio_jaintianensis',
                'Danio_tinwini']:
    scbuscos[species] = set(os.listdir(f'../../data/busco-out/{species}/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences'))

scbuscoset = set.intersection(*scbuscos.values())

for species in ['Danio_albolineatus',
                'Danio_choprai',
                'Danio_jaintianensis',
                'Danio_tinwini']:
    with open(f'../../data/danio-reads/{species}_busco.bed', 'w') as outfile:
        for bid in list(scbuscoset)[:20]:
            record = list(SeqIO.parse(f'../../data/busco-out/{species}/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/{bid}',
                'fasta'))[0]
            chrom, start, stop = re.split(':|-', record.id)
            if int(start) > int(stop):
                start, stop = stop, start
            line = [chrom, start, stop, species + '_' + bid.strip('.faa')]
            outfile.write('\t'.join(line) + '\n')


