#!/usr/bin/env python3

import sys

data = []
with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith('#'):
            data.append(line.strip())
            continue
        line = line.strip().split('\t')
        chrom = line[0].split(':')[0]
        blockstart = int(line[0].split(':')[1].split('-')[0])
        line[0] = chrom
        line[3] = str(int(line[3]) + blockstart)
        line[4] = str(int(line[4]) + blockstart)
        line = '\t'.join(line)
        data.append(line)

with open(sys.argv[2], 'w') as outfile:
    for line in data:
        outfile.write(f'{line}\n')
        

    
