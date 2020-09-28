#!/usr/bin/env python3

import argparse
import re
import numpy as np
import os

desc='Generate pairwise distance matrix in phylip format'
parser=argparse.ArgumentParser(description=desc)
parser.add_argument('input')
parser.add_argument('output')
args = parser.parse_args()

def kimura_pdist(m, npos, gaps):
    D = uncorrected_dist(m, npos, gaps)
    return -np.log(1 - D - (0.2*(D**2)))

def uncorrected_dist(m, npos, gaps):
    return 1-m/(npos-gaps)

def relabel(label):
    parts = label.split('.')[0].split('_')
    gcode = parts[0][0]
    scode = ''.join(parts[1][:3])
    gene = parts[2]
    return f'{gcode}{scode}_{gene}'

def parse_alignments(filename):
    pattern = ('# 1: (.+)\n', 
               '# 2: (.+)\n',
               '.+\n.+\n.+\n.+\n',
               '# Length: (\d+)\n',
               '# Identity:\s+(\d+).+\n',
               '.+\n',
               '# Gaps:\s+(\d+)')
    pattern = re.compile(''.join(pattern), re.MULTILINE)
    distance = {}
    
    with open(filename) as infile:
        data = infile.read()
        hits = re.findall(pattern, data)
        for align in hits:
            g1, g2, npos, m, gaps = align
            g1, g2 = relabel(g1), relabel(g2)
            distance[(g1, g2)] = kimura_pdist(int(m), int(npos), int(gaps))
            # distance[(g1, g2)] = uncorrected_dist(int(m), int(npos), int(gaps))
    return distance    
    # if set(i[0] for i in distance.keys()) == set(i[1] for i in distance.keys()):
    #     return distance
    # else:
    #     raise ValueError('incorrect distance matrix')

def main():
    distance = {}
    c = 0
    for filename in os.listdir(args.input):
        if re.match('needle\d+.out', filename):
            c += 1
            if c % 10 == 0:
                print(c, '...')
            distance.update(parse_alignments(f'{args.input}/{filename}'))
    geneidsi = sorted(set(i[0] for i in distance.keys()))
    geneidsj = sorted(set(i[1] for i in distance.keys()))
    print(len(geneidsi), len(geneidsj)) 
    with open(args.output, 'w') as outfile:
        datlen = str(len(geneidsi))
        data = f'{datlen:>5}\n'
        for i in geneidsi:
            line = f'{i:10}  '
            for j in geneidsj:
                # if j >= i:
                #     continue
                average = np.mean([distance[(i, j)], distance[(j, i)]])
                line += f'\t{average:.3f}'
            line = line + '\n'
            data += line
        data = data
        outfile.write(data)

if __name__ == '__main__':
    main()
        
