#!/usr/bin/env python3

from collections import defaultdict
import re
import sys

def load_clusters(filename):
    clusters = defaultdict(list)
    with open(filename) as infile:
        for line in infile:
            if line.startswith('>'):
                key = re.match('>Cluster (\d+)', line)
            else:
                seqid = re.search('[\w\s]+, >(.+)\.{3}', line).group(1)
                clusters[key.group(1)].append(seqid)
    return clusters

def parse_clusters(clusters):
    cluster_types = defaultdict(list) 
    for key, val in clusters.items():
        if len(val) == 1:
            cluster_types['singleton'].append(key)
            continue
        splist, spset = [], set()
        for seqid in val:
            species = re.match('(\w+_\w+)_', seqid).group(1)
            splist.append(species)
            spset.add(species)
        if len(spset) == len(splist):
            cluster_types['orthologs'].append(key)
        elif len(spset) == 1:
            cluster_types['paralogs'].append(key)
        else:
            cluster_types['mixed'].append(key)
    return cluster_types


if __name__ == '__main__':
    c = load_clusters(sys.argv[1])
    ct = parse_clusters(c)
    for key in ct[sys.argv[2]][::-1]:
        print(key)
        for seqid in c[key]:
            print(f'\t{seqid}')


