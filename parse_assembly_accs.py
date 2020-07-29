#!/usr/bin/env python3

import xml.etree.ElementTree as et

tree = et.parse('/Users/jonwells/Genomes/Cypriniformes/assembly_result.xml')
root = tree.getroot()

for docsum in root:
    acc = docsum.find('AssemblyAccession')
    org = docsum.find('Organism')
    print(f'{acc.text}\t{org.text}')
