#!/usr/bin/env python3

import re
import pandas as pd

##### LEGACY #####


# This script checks overlaps of denovo finz annotations with those of ensembl,
# and renames denovo with existing Ensembl annotations, where appropriate.
# Relies on the output of gffcompare -r ensembl_finz_znf.gff denovo_finz_znf.gff

gffcmp_df = pd.read_csv('../../data/gffs/gffcmp.denovo_finz_znf.gff.tmap', sep='\t')
keep_df = gffcmp_df.loc[~gffcmp_df.class_code.isin(['c', '='])]
rename_df = gffcmp_df.loc[gffcmp_df.class_code.isin(['c', '='])]

# Use these to rename with ensembl ids
# transcript_dict = dict(zip(rename_df['qry_id'], rename_df['ref_id']))
# gene_dict = dict(zip(rename_df['qry_gene_id'], rename_df['ref_gene_id']))

# Use these to keep names the same
transcript_dict = dict(zip(rename_df['qry_id'], rename_df['qry_id']))
gene_dict = dict(zip(rename_df['qry_gene_id'], rename_df['qry_gene_id']))

with open('../../data/gffs/Danio_rerio_finz.final.gff', 'w') as outfile:
    with open('../../data/gffs/denovo_finz_znf.gff') as infile:
        for line in infile:
            line = line.strip().split('\t')
            if line[2] == 'gene':
                gid = re.match('ID=(\w+)$', line[-1]).group(1)
                id_replacement = gene_dict.get(gid, gid)
                line[-1] = f'ID={id_replacement}'
            elif line[2] == 'transcript':
                gid = re.match('ID=(\w+\.t1);Parent=(\w+)$', line[-1])
                id_replacement = transcript_dict.get(gid.group(1), gid.group(1)).strip('transcript:')
                parent_replacement = gene_dict.get(gid.group(2), gid.group(2))
                line[-1] = f'ID={id_replacement};Parent={parent_replacement}'
                line[2] = 'mRNA'
            elif line[2] in ['start_codon', 'stop_codon']:
                gid = re.match('Parent=(\w+\.t1)$', line[-1]).group(1)
                parent_replacement = transcript_dict.get(gid, gid).strip('transcript:')
                line[-1] = f'Parent={parent_replacement}'
            elif line[2] == 'CDS':
                gid = re.match('ID=(\w+\.t1)\.cds;Parent=(\w+\.t1)$', line[-1])
                id_replacement = transcript_dict.get(gid.group(1), gid.group(1)).strip('transcript:')
                parent_replacement = transcript_dict.get(gid.group(2), gid.group(2)).strip('transcript:')
                line[2] = 'exon'
                line[-1] = f'ID={id_replacement}.exon;Parent={parent_replacement}'
                outfile.write('\t'.join(line) + '\n')
                line[2] = 'CDS'
                line[-1] = f'ID={id_replacement}.cds;Parent={parent_replacement}'
            outfile.write('\t'.join(line) + '\n')

## Best to remove all the "transcript:" and "gene:" after with sed
# sed 's/transcript://g'  Danio_rerio.GRCz11.101.curated_finz.gff3 | sed 's/gene://g'
