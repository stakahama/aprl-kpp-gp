#!/usr/bin/env python

################################################################################
## runs executables
##
## ~exec_dual.py~
##
##
## satoshi.takahama@epfl.ch
##
################################################################################


###_* -------------------- import libraries --------------------

import os
import subprocess
from collections import OrderedDict
import pandas as pd
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='generate a0.txt')
parser.add_argument('ROOT',type=str)
parser.add_argument('NUMBERS',type=str)
args = dict(vars(parser.parse_args()))

###_* -------------------- define paths --------------------

args['MAINPATH'] = os.path.dirname(__file__)
HERE = os.getcwd()
paths = OrderedDict([(p, os.path.join(HERE,'exec_'+p))  for p in ['gas','total']])

###_* -------------------- parse numbers --------------------

if args['NUMBERS']=='all':
    pass
else:
    numbers = map(int,args['NUMBERS'].split(','))

###_* -------------------- read output --------------------

vp = pd.read_csv('{ROOT}_props_298.csv'.format(**args),index_col='compound')

for i in numbers:
    runpath = 'run_{:03d}'.format(i)
    output = pd.read_csv(os.path.join(HERE,runpath,'gas','{ROOT}_formatted.txt'.format(**args)),
                         index_col='TIME')
    
