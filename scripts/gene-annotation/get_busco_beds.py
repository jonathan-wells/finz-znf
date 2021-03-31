#!/usr/bin/env python

import os

scbuscos = {}
for species in ['Danio_albolineatus',
                'Danio_choprai',
                'Danio_jaintianensis',
                'Danio_tinwini']:
    scbuscos[species] = os.listdir(f'../../data/busco-out/{species}/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences')

print(scbuscos)
