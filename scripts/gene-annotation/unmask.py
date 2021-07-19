#!/usr/bin/env python3

from Bio import SeqIO, Seq
import sys

def unmask(fastafile, blastfile, outfilename):
    records = SeqIO.to_dict(SeqIO.parse(fastafile, 'fasta'))
    with open(blastfile) as infile:
        for line in infile:
            line = line.strip().split('\t')
            chrom, sstart, send = line[1], int(line[8]), int(line[9])
            sstart, send = min([sstart, send]) - 1, max([sstart, send])
            currseq = str(records[chrom].seq)
            currseq = currseq[:sstart] + currseq[sstart:send].upper() + currseq[send:]
            records[chrom].seq = Seq.Seq(currseq)
    SeqIO.write(records.values(), outfilename, 'fasta')

if __name__ == '__main__':
    unmask(sys.argv[1], sys.argv[2], sys.argv[3])
