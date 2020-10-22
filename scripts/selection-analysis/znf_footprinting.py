#!/usr/bin/env python3

import re
from Bio import SeqIO
from collections import defaultdict

def footprint(seq):
    znf1 = re.compile('C.{2,4}C.{8,15}H.{3,5}H')
    fp_index = (-12, -10, -9, -6)
    hits = re.findall(znf1, seq)
    
    footprint = []
    for hit in hits:
        footprint.append(''.join([hit[i] for i in fp_index]))
    return tuple(footprint)

def extract_footprints(min_znfs=3):
    footprints = {}
    for record in SeqIO.parse('../../data/seqs/cypriniformes_augustus_finz.fa',
                              'fasta'):
        fp = footprint(str(record.seq))
        if len(fp) >= min_znfs:
            footprints[record.id] = fp
    return footprints

if __name__ == '__main__':
    fps = extract_footprints()
    all_fps = defaultdict(int)
    for key, val in fps.items():
        pretty_fp = ','.join(val)
        # print(f'>{key}')
        # print(pretty_fp)
        # print(f'{key}\t{pretty_fp}')
        for i in val:
            all_fps[i] += 1

    for i in sorted(all_fps.items(), key=lambda x: x[1]):
        print(i[0], i[1])
    print(len(all_fps))
    print(sum(all_fps.values()))

    # fps_count = defaultdict(int)
    # for val in fps.values():
    #     fps_count[val] += 1
    # for i in sorted(fps_count.items(), key=lambda x: x[1]):
    #     print(i[0], i[1])
