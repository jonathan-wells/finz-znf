#!/usr/bin/env python3

from collections import defaultdict
import re
import sys
from Bio import SeqIO
import argparse

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

def get_seqs(clusters, cid, codons=False):
    if codons:
        seqs = SeqIO.to_dict(SeqIO.parse('../../data/seqs/cypriniformes_augustus_finz.dna.fa', 'fasta'))
    else:
        seqs = SeqIO.to_dict(SeqIO.parse('../../data/seqs/cypriniformes_augustus_finz.fa', 'fasta'))
    for i in clusters[str(cid)]:
        SeqIO.write(seqs[i], sys.stdout, 'fasta')
    

def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(filename='../../data/cdhit_clusters.clstr')
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       dest='command',
                                       help='additional help')
    
    parse = subparsers.add_parser('parse')
    parse.add_argument('-f', '--filename', help='cd-hit output file', type=str)
    parse.add_argument('-t', '--type',
            required=True,
            type=str,
            help='Cluster type: choose from "orthologs", "paralogs", "mixed"')
    parse.set_defaults(func=parse_clusters)
    
    getseqs = subparsers.add_parser('getseqs')
    getseqs.add_argument('-f', '--filename', help='cd-hit output file', type=str)
    getseqs.add_argument('-c', '--codons', type=bool)
    getseqs.add_argument('-i', '--id', type=str)
    getseqs.set_defaults(func=get_seqs)
    
    args = parser.parse_args()
    clusters = load_clusters(args.filename)
    if args.command == 'parse':
        parsed_clusters = args.func(clusters)
        for key in parsed_clusters[args.type][::-1]:
            print(key)
            for seqid in clusters[key]:
                print(f'\t{seqid}')
    elif args.command == 'getseqs':
        if args.codons:
            args.func(clusters, args.id, args.codons)
        else:
            args.func(clusters, args.id)
        


if __name__ == '__main__':
    main()
