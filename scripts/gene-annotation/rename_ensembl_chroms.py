#!/usr/bin/env python3

import sys

chr2acc = {}
with open('../../data/misc/GCF_000002035.6_GRCz11_assembly_report.txt') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.split()
        if line[3] == 'Chromosome' or line[3] == 'Mitochondrion':
            chr2acc[line[2]] = line[6]
        elif line[3] == 'na':
            chr2acc[line[4]] = line[6]


data = []
with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith('#'):
            data.append(line.strip())
            continue
        line = line.strip().split('\t')
        line[0] = chr2acc[line[0]]
        line = '\t'.join(line)
        data.append(line)

with open(sys.argv[2], 'w') as outfile:
    for line in data:
        outfile.write(f'{line}\n')
