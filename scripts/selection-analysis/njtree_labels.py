#!/usr/bin/env python3

""" Because of the peculiar requirements of Phylip's distance matrix format, the
    resulting labels in the Neighbor-Joining tree are not very informative. This
    short script corrects this issue.
"""

from Bio import SeqIO

def relabel(label):
    parts = label.split('.')[0].split('_')
    gcode = parts[0][0]
    scode = ''.join(parts[1][:3])
    gene = parts[2]
    return f'{gcode}{scode}_{gene}'

labeldict = {}
for record in SeqIO.parse('../../data/seqs/cypriniformes_augustus_finz.fa', 'fasta'):
    njtree_label = relabel(record.id) + ':'
    newlabel = f'{record.id}_[len={len(record.seq)}]:'
    labeldict[njtree_label] = newlabel

with open('../../data/selection-analysis/needle_nj_incgaps_relabeled.nwk', 'w') as outfile:
    with open('../../data/selection-analysis/needle_nj_incgaps.nwk') as infile:
        newtree = infile.read().replace('\n', '')
        for key, val in labeldict.items():
            newtree = newtree.replace(key, val)
        outfile.write(newtree)

