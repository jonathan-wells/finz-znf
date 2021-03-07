#!/usr/bin/env python3

import sys
import re
from collections import defaultdict

# Convert from .bed to .gtf format suitable for TEtranscripts. Important
# features include omission of all low_complexity type repeats, and naming of
# individual transcripts with dup1..n

# GTF format:
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

def ucsc2ncbi():
    accdict = {}
    with open('../../data/misc/GCF_000002035.6_GRCz11_assembly_report.txt') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            accdict[line[9]] = line[6]
    return accdict

def parse_attrs(line, dupdict):
    gene_id = line[3]
    dupnum = dupdict.get(gene_id, 0)
    if dupnum == 0:
        transcript_id = line[3]
    else:
        transcript_id = f'{line[3]}_dup{dupnum}'
    clsfam = line[10].split('/')
    if len(clsfam) == 1:
        family_id, class_id = clsfam[0], clsfam[0]
    elif len(clsfam) == 2:
        class_id, family_id = clsfam
    else:
        raise AttributeError(f'malformatted class: {line[10]}')
    return transcript_id, gene_id, family_id, class_id

def format_attrs(transcript_id, gene_id, family_id, class_id):
    return f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; family_id "{family_id}"; class_id "{class_id}";'

def main():
    dupdict = defaultdict(int)
    accdict = ucsc2ncbi()
    with open(sys.argv[2], 'w') as outfile:
        with open(sys.argv[1]) as infile:
            for line in infile:
                line = line.strip().split('\t')
                seqname = accdict[line[0]]
                source = 'RepeatMasker'
                feature = 'exon'
                start = int(line[1]) + 1  # GFF/GTF uses 1-based indexing
                end = int(line[2]) + 1 
                score = line[4]
                strand = line[5]
                frame = '.'
                transcript_id, gene_id, family_id, class_id = parse_attrs(line, dupdict)

                if class_id in ['Low_complexity', 'Retroposon', 'Simple_repeat', 'rRNA', 'scRNA', 'snRNA', 'tRNA']:
                    continue
                elif gene_id in ['SAT-7_DR', 'BEL-64_DRe-I']:
                    # These are not TEs but FINZ-ZNFs, ironically.
                    continue
                attrs = format_attrs(transcript_id, gene_id, family_id, class_id)
                newline = f'{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attrs}\n'
                outfile.write(newline)
                dupdict[gene_id] += 1

if __name__ == '__main__':
    main()

