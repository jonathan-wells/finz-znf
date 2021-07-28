#!/usr/bin/env python

from collections import defaultdict
import re
from Bio import SeqIO, Seq


class Annotate:

    def __init__(self, fastafile, pamlfile=None, rstfile=None):
        self.records = list(SeqIO.parse(fastafile, 'fasta'))
        self.c2h2_pattern = r'(C)..(C).{2,15}(..).(..).(..)(H)...(H)'
        self.psel_pattern = r'\s*(\d+)\s+[A-Z]\s+(0.\d+)\**\s+(\d+\.\d+)'
        self.rst_neb_pattern = r'\s*(\d+)\s+[A-Z].+\(\s*[1-9]+\)\s+(\d+\.\d+)\s+(\d+\.\d+)'
        self.rst_beb_pattern = r'\s*(\d+)\s+[A-Z].+(\d+\.\d+)\s+\(\s*[1-9]+\)\s+(\d+\.\d+)\s+\+\-'
        if pamlfile:
            self.paml = self._parse_paml(pamlfile) 
        if rstfile:
            self.rst = self._parse_rst(rstfile)

    def annotate_c2h2(self, outfilename, fmt='jff'):
        with open(outfilename, 'w') as outfile:
            if fmt == 'jff':
                data = self._annotate_c2h2_jff()
            elif fmt == 'gff':
                data = self._annotate_c2h2_gff()
            for line in data:
                outfile.write(f'{line}\n')
    
    def annotate_positively_selected(self, outfilename, model, fmt='jff'):
        with open(outfilename, 'w') as outfile:
            if fmt == 'jff':
                data = self._annotate_psel_jff(model)
            elif fmt == 'gff':
                data = self._annotate_psel_gff(model)
            for line in data:
                outfile.write(f'{line}\n')
    
    def annotate_rst(self, outfilename, model, fmt='gff'):
        with open(outfilename, 'w') as outfile:
            if fmt == 'jff':
                raise KeyError('jff format not yet implemented for rst output')
            elif fmt == 'gff':
                data = self._annotate_rst_gff(model)
            for line in data:
                outfile.write(f'{line}\n')

    def _parse_paml(self, pamlfile):
        pamldict = defaultdict(list)
        with open(pamlfile) as infile:
            for line in infile:
                if re.match('^Naive Empirical Bayes \(NEB\) analysis$', line):
                    nebflag = True
                elif re.match('^Bayes Empirical Bayes \(BEB\) analysis \(Yang', line):
                    nebflag = False
                elif re.match('^Model 2:', line):
                    m12flag = True
                elif re.match('^Model 8:', line):
                    m12flag = False
        
                psel_sites = re.search(self.psel_pattern, line)
                if not psel_sites:
                    continue
                if nebflag and m12flag:
                    pamldict['neb12'].append(psel_sites.groups())
                elif nebflag and not m12flag:
                    pamldict['neb78'].append(psel_sites.groups())
                elif not nebflag and m12flag:
                    pamldict['beb12'].append(psel_sites.groups())
                elif not nebflag and not m12flag:
                    pamldict['beb78'].append(psel_sites.groups())
        return pamldict

    def _parse_rst(self, rstfile):
        rstdict = defaultdict(list)
        with open(rstfile) as infile:
            for line in infile:
                if re.match('^Naive Empirical Bayes \(NEB\) probabilities', line):
                    nebflag = True
                elif re.match('^Bayes Empirical Bayes \(BEB\) probabilities', line):
                    nebflag = False
                elif re.match('^Model 2:', line):
                    m12flag = True
                elif re.match('^Model 8:', line):
                    m12flag = False
        
                rst_neb_sites = re.search(self.rst_neb_pattern, line)
                rst_beb_sites = re.search(self.rst_beb_pattern, line)
                if not (rst_neb_sites or rst_beb_sites):
                    continue
                if nebflag and m12flag:
                    rstdict['neb12'].append(rst_neb_sites.groups())
                elif nebflag and not m12flag:
                    rstdict['neb78'].append(rst_neb_sites.groups())
                elif not nebflag and m12flag:
                    bsites = rst_beb_sites.groups()
                    rstdict['beb12'].append((bsites[0], bsites[2], bsites[1]))
                elif not nebflag and not m12flag:
                    bsites = rst_beb_sites.groups()
                    rstdict['beb78'].append((bsites[0], bsites[2], bsites[1]))
        return rstdict
                    
    def _annotate_psel_jff(self, model):
        data = []
        for record in self.records:
            for site in self.paml[model]:
                data.append(f'seq\t{record.id}\t-1\t{site[0]}\t{site[0]}\tPsel\t{site[1]}')
        return data

    def _annotate_psel_gff(self, model):
        data = ['##gff-version 3']
        psel_sites = [int(i[0]) for i in self.paml[model]]
        for record in self.records[0:1]:
            for i in range(1,len(record.seq.translate())+1):
                if i in psel_sites:
                    continue
                data.append(f'{record.id}\tseq\tbase\t{i}\t{i}\t0.0\t+\t.\t.')
            for site in self.paml[model]:
                data.append(f'{record.id}\tseq\tPsel\t{site[0]}\t{site[0]}\t{site[1]}\t+\t.\t.')
        data = ['##gff-version 3'] + list(set(data))
        return data

    def _annotate_c2h2_jff(self):
        data = []
        for record in self.records:
            translation = str(record.seq.translate())
            locs = re.finditer(self.c2h2_pattern, translation)
            for m in locs:
                data.append(f'seq\t{record.id}\t-1\t{m.start(1)+1}\t{m.end(1)}\tCys\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(2)+1}\t{m.end(2)}\tCys\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(3)+1}\t{m.end(3)}\tBase-contacting\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(4)+1}\t{m.end(4)}\tBase-contacting\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(5)+1}\t{m.end(5)}\tBase-contacting\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(6)+1}\t{m.end(6)}\tHis\t0.0')
                data.append(f'seq\t{record.id}\t-1\t{m.start(7)+1}\t{m.end(7)}\tHis\t0.0')
        return data
    
    def _annotate_c2h2_gff(self):
        data = ['##gff-version 3']
        for record in self.records[0:1]:
            translation = str(record.seq.translate())

            locs = list(re.finditer(self.c2h2_pattern, translation))
            for m in locs:
                data.append(f'{record.id}\tseq\tCys\t{m.start(1)+1}\t{m.end(1)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tCys\t{m.start(2)+1}\t{m.end(2)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(3)+1}\t{m.end(3)-1}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(3)+2}\t{m.end(3)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(4)+1}\t{m.end(4)-1}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(4)+2}\t{m.end(4)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(5)+1}\t{m.end(5)-1}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tBase-contacting\t{m.start(5)+2}\t{m.end(5)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tHis\t{m.start(6)+1}\t{m.end(6)}\t0.0\t+\t.\t.')
                data.append(f'{record.id}\tseq\tHis\t{m.start(7)+1}\t{m.end(7)}\t0.0\t+\t.\t.')
            
            c2h2_sites = [m.start(i) + 1 for i in range(1, 8) for m in locs] + \
                         [m.start(i) + 2 for i in range(3, 6) for m in locs]
            for i in range(1,len(translation)+1):
                if i in c2h2_sites:
                    continue
                data.append(f'{record.id}\tseq\tbase\t{i}\t{i}\t0.0\t+\t.\t.')
        return data

    def _annotate_rst_gff(self, model):
        data = ['##gff-version 3']
        record = self.records[0]
        for site in self.rst[model]:
            data.append(f'{record.id}\tseq\tPsel\t{site[0]}\t{site[0]}\t{site[1]}\t+\t.\tPr:{site[2]}')
        return data


def main():
    for i in range(1,8):
        cdir = f'../../data/selection-analysis/clade{i}'
        clade = Annotate(f'{cdir}/pairwise_needleman_clade{i}.final.fa',
                         f'{cdir}/pairwise_needleman_clade{i}.final.out',
                         f'{cdir}/rst')
        clade.annotate_c2h2(f'{cdir}/clade{i}_c2h2.jff', 'jff')
        clade.annotate_positively_selected(f'{cdir}/clade{i}_psel_neb12.jff', 'neb12', 'jff')
        clade.annotate_c2h2(f'{cdir}/clade{i}_c2h2.gff', 'gff')
        clade.annotate_positively_selected(f'{cdir}/clade{i}_psel_neb12.gff', 'neb12', 'gff')
        clade.annotate_rst(f'{cdir}/clade{i}_rst_neb12.gff', 'neb12', 'gff')
    
if __name__ == '__main__':
    main()
