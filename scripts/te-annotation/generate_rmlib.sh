#!/usr/bin/env bash

# Have to use older repbase library until someone sorts out classification
# format discrepancies.
/programs/RepeatMasker/util/buildRMLibFromEMBL.pl \
    /programs/RepeatMasker_4-0-8/Libraries/RepeatMaskerLib.embl \
    > ../../data/rmlibs/rmlib.fa

