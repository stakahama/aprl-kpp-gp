#!/usr/bin/env python

################################################################################
## runs executables
##
## ~reinitialize_a0.py~
##
## 
##
## satoshi.takahama@epfl.ch
##
################################################################################


###_* -------------------- import libraries --------------------

import os
import sys
scriptspath = os.path.join(os.path.dirname(os.path.dirname(__file__)),'scripts')
sys.path.append(scriptspath)
from simpol2 import SimpolClass
from collections import OrderedDict
import subprocess
import pandas as pd
import numpy as np
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='generate a0.txt')
parser.add_argument('ROOT',type=str)
parser.add_argument('INPUTRUN',type=str)
parser.add_argument('OUTPUTRUN',type=str)
args = dict(vars(parser.parse_args()))

###_* -------------------- define paths --------------------

a0filename = 'molefrac_init.txt'

args['MAINPATH'] = os.path.dirname(__file__)
HERE = os.getcwd()
paths = OrderedDict([(p, os.path.join(HERE,'exec_'+p))  for p in ['gas','total']])

###_* -------------------- calculate partitioning --------------------

###_ . define partitioning functions and vapor pressures

###_ . define constant and partitioning function
cfactor = 10**9          # atm to ppb conversion factor

def partition(p,p0):     # p0, p in ppb
    excess = p-p0
    aer = excess.map(lambda x: x if x > 0 else 0)
    return aer/aer.sum() # mole fraction

###_ . read in temperature
with open('{ROOT}.def'.format(**args)) as f:
    for line in f:
        if 'TEMP' in line:
            temp = float(line.split('=')[1].strip())

###_ . set in/out runpaths
runpath = OrderedDict()
runpath['INP'] = os.path.join(HERE,args['INPUTRUN'])
runpath['OUT'] = os.path.join(HERE,args['OUTPUTRUN'])

if not os.path.exists(runpath['OUT']):
    os.mkdir(runpath['OUT'])

resultsfile = os.path.join(runpath['INP'],'gas/{ROOT}_formatted.csv'.format(**args))
a0file = os.path.join(runpath['OUT'],a0filename)

###_ . calculate vapor pressures
simp = SimpolClass()
simp.read_compounds('{ROOT}_SIMPOLGroups.csv'.format(**args))
vp = simp.calc_properties(temp)
# vp = pd.read_csv('{ROOT}_props_298.csv'.format(**args),index_col='compound')

###_ . calculate mole fractions
last = pd.read_csv(resultsfile,index_col='TIME').iloc[-1]
molefrac = pd.DataFrame(partition(last.ix[vp.index],vp['p0']*cfactor),
                        columns=['molefrac'])

###_ . merge with index
indices = pd.read_csv('compound_indices_table.csv',index_col='compound')
table = indices[['index']].join(molefrac,how='inner').sort('index')

###_ . export
table.to_csv(a0file,sep='\t',header=False,index=False)
