#!/usr/bin/env python3

import sys

with open('../../data/misc/chr2acc.txt') as infile:
    infile.readline()
    chr2acc = {i.split()[0]: i.split()[1] for i in infile}

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

with open(sys.argv[1], 'w') as outfile:
    for line in data:
        outfile.write(f'{line}\n')
