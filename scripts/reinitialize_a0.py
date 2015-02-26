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
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='run kpp for gas-phase and total simulation')
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

{'{ROOT}_formatted.csv'
