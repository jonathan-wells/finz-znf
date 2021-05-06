#!/usr/bin/env python

import treeswift as ts
import re

def load_tree(filename):
    tf = open(filename) \
            .read() \
            .replace('\n', '')
    tree = ts.read_tree(tf, 'Newick')
    for node in tree.traverse_postorder():
        if node.is_root():
            root = node
            return tree, root

def length_mismatch(clade):
    lengths = []
    for label in clade:
        lengths.append(int(re.search('len=(\d+)', label).group(1)))
    return max(lengths) - min(lengths)

def get_clades(tree, root, minleaves=10,maxmismatch=100,species='Danio_rerio'):
    clades = []
    for node in root.traverse_bfs():
        leaves = [leaf.get_label() for leaf in node[0].traverse_leaves()]
        if len(leaves) < minleaves:
            continue
        paralogues = set(leaf for leaf in leaves if leaf.startswith(species))
        if len(paralogues) == len(leaves) and length_mismatch(paralogues) <= maxmismatch:
            clades.append(paralogues)
    
    # Remove clades that are subsets of others in list
    subsets = set()
    clades = {i: clade for i, clade in enumerate(sorted(clades, key=len))}
    for key1, val1 in clades.items():
        for key2, val2 in clades.items():
            if key1 == key2:
                continue
            if val1.issubset(val2):
                subsets.add(key1)
    for key in list(subsets):
        clades.pop(key)
    return clades

if __name__ == '__main__':
    tree, root = load_tree('../../data/selection-analysis/pairwise_needleman_longlabel.nwk')
    clades = get_clades(tree, root, 10, 30, 'Danio_rerio')
    i = 1
    for val in clades.values():
        with open(f'../../data/selection-analysis/pairwise_needleman_clade{i}.taxa', 'w') as outfile:
            outfile.write('\n'.join(sorted(val)))
        i += 1

    
