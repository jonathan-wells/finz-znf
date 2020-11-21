#!/usr/bin/env python3

import sys

chr2acc = {}
with open('../../data/misc/GCF_000002035.6_GRCz11_assembly_report.txt') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.split()
        if line[3] == 'Chromosome':
            chr2acc[line[-1].strip('chr')] = line[6]
        elif line[3] == 'Mitochondrion':
            chr2acc[line[2]] = line[6]
        elif line[3] == 'na':
            chr2acc[line[4]] = line[6]


data = []
with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith('##gff-version 3'):
            data.append(line.strip())
        elif line.startswith('##sequence-region'):
            line = line.split()
            line[1] = chr2acc[line[1]]
            line = f'{line[0]}   {line[1]} {line[2]} {line[3]}'
            data.append(line)
        elif line.startswith('#!'):
            data.append(line.strip())
        elif line.startswith('###'):
            continue
        else:
            line = line.strip().split('\t')
            line[0] = chr2acc[line[0]]
            line = '\t'.join(line)
            data.append(line)

with open(sys.argv[2], 'w') as outfile:
    for line in data:
        outfile.write(f'{line}\n')
