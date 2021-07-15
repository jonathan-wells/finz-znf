#!/usr/bin/env python3

"""This script takes AUGUSTUS gff files as input, and offsets the genomic
according to requirements - relative either to the start of chromosomes, or the 
start of genes."""

import sys
import re

def offset_gff(filename):
    """Return gff offset to true chromosome starts"""
    data = []
    with open(filename) as infile:
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
    return data

def offset_gff2(filename, flank=0, prefix=None):
    """Return gff offset to length of gene."""
    data = []
    with open(filename) as infile:
        for line in infile:
            if line.startswith('#'):
                data.append(line.strip())
                continue
            line = line.strip().split('\t')
            if line[2] == 'gene':
                offset = int(line[3])
                chrom = prefix + re.match('ID=(\w+)', line[-1]).group(1)
            line[0] =  chrom
            line[3] = str(int(line[3]) - offset + flank + 1)
            line[4] = str(int(line[4]) - offset + flank + 1)
            data.append('\t'.join(line))
    return data

def offset_gff3(filename, blockfilename):
    data = []
    blocks = set()
    with open(blockfilename) as blockfile:
        for blockline in blockfile:
            blockline = blockline.strip().split('\t')
            blockstart, blockend = int(blockline[1]), int(blockline[2])
            blockchrom = f'{blockline[0]}:{blockstart}-{blockend}'
            blocks.add((blockchrom, blockstart, blockend))
    with open(filename) as infile:
        for line in infile:
            if line.startswith('#'):
                data.append(line.strip())
                continue
            line = line.strip().split('\t')
            for block in blocks:
                blockchrom, blockstart, blockend = block
                if (int(line[3]) > blockstart and int(line[3]) < blockend) and (int(line[4]) > blockstart and int(line[4]) < blockend):
                    print(blockchrom)
                    line[0] = blockchrom
                    line[3] = str(int(line[3]) - blockstart)
                    line[4] = str(int(line[4]) - blockstart)
                    break
            data.append('\t'.join(line))
    return data
                    

if __name__ == '__main__':
    if len(sys.argv) == 3:
        data = offset_gff(sys.argv[1])    
        with open(sys.argv[2], 'w') as outfile:
            for line in data:
                outfile.write(f'{line}\n')
    elif len(sys.argv) == 5:
        data = offset_gff2(sys.argv[1], int(sys.argv[3]), sys.argv[4])    
        with open(sys.argv[2], 'w') as outfile:
            for line in data:
                outfile.write(f'{line}\n')
    elif len(sys.argv) == 4:
        print('fixing_hints')
        data = offset_gff3(sys.argv[1], sys.argv[2])
        with open(sys.argv[3], 'w') as outfile:
            for line in data:
                outfile.write(f'{line}\n')
    else:
        raise ValueError('Wrong number of arguments, check script usage.')
    
