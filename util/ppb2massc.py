#!/usr/bin/env python

################################################################################
## convert units
##
## ~ppb2massc.py~
##
## USAGE:
## $ python [...]
##
## EXAMPLES:
## $ python [...]
##
## satoshi.takahama@epfl.ch
##
################################################################################

###_* -------------------- import libraries --------------------

import os
import numpy as np
import pandas as pd
from functools import partial
from operator import mul
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='convert _formatted.csv (ppb) to _org_mass.csv (micrograms per cubic meter)')
parser.add_argument('ROOT',type=str)
parser.add_argument('RUNPATH',type=str)
args = dict(vars(parser.parse_args()))

def ppb2massc(ppb, molwt, temp=298.15, press=1):
    R = 8.206e-5
    return 1e-3*ppb*molwt/R*press/temp

###_* -------------------- read molecular weights --------------------

molwt = pd.read_table('mcm_{ROOT}_mass.txt'.format(**args),
                      skiprows=18,header=None,skipinitialspace=True,
                      names=['compound','SMILES','InChI','molwt'],
                      index_col='compound').molwt
molwt.index = molwt.index.map(str.strip)

###_* -------------------- read/convert ppb files --------------------

folders = ['gas','total']
inpfiles = [
    '{ROOT}_formatted.csv'.format(**args),
    '{ROOT}_aer_formatted.csv'.format(**args)
    ]
newfile = lambda x: x.replace('_formatted.csv','_org_massc.csv')

for d in folders:
    for f in inpfiles:
        ## create file name
        filename = os.path.join(args['RUNPATH'],d,f)
        if not os.path.exists(filename):
            continue
        ## import and convert
        ppb = pd.read_csv(filename,index_col='TIME')
        massc = ppb[molwt.index].apply(partial(ppb2massc,molwt=molwt),axis=1)
        ## export
        massc.reset_index().to_csv(newfile(filename),index=False)
        print '...', newfile(filename), 'created'
