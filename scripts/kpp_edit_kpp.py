# -*- coding: utf-8 -*-
"""
This script
1) adds stoichiometric equations to orgainc compounds
2) fixes equation syntax (exponential equations)
3) removes inorganic compounds (they are handled through inorganic.kpp)
4) dos2unix

output of this file and editkpp.py+dos2unix are identical according to 'diff -q'.

* this script works but should be rewritten (low priority)
- SmilesMatch can reuse code from structmatch


Created on Mon Aug 11 16:49:52 2014
Edited Feb 2015
@author: shipley, stakahama
"""

###_* -------------------- import libraries --------------------

import os
import pandas as pd
import numpy as np
import pybel
from collections import OrderedDict
from operator import itemgetter
from kpp_generate_SIMPOLGroups import kppParameters
import argparse

###_* -------------------- parse arguments --------------------

# root = raw_input("ROOT: ")
parser = argparse.ArgumentParser(description='edit {root}.kpp file')
parser.add_argument('root',type=str)
args = parser.parse_args()
root = args.root

###_* -------------------- read mass_{ROOT}_mass.txt --------------------

class SmilesMatch:

    def __init__(self):
        ## SMARTS patterns required by SIMPOL
        self.smartspatt = OrderedDict([
            ('C','[#6]'),
            ('O','[#8]'),
            ('H','nomatch'),
            ('N','[#7]'),
            ('S','[#16]'),
            ])
        self.smartsaux = OrderedDict([
            ('h1','[*h1]'),
            ('h2','[*h2]'),
            ('h3','[*h3]'),
            ('h4','[*h4]'),
            ])

    def get_abundances(self,smilesstr=None):
        ## main body
        mol = pybel.readstring('smi',smilesstr)
        ## store SIMPOL patterns in ordered dictionary
        abundances = OrderedDict()
        for key,patt in self.smartspatt.items()[0:]:
            abundances[key] = 0 if 'nomatch' in patt else \
                              len(pybel.Smarts(patt).findall(mol))
        ## find auxiliary patterns
        aux = OrderedDict()
        for key,patt in self.smartsaux.items():
            aux[key] = len(pybel.Smarts(patt).findall(mol))
        ## combine and return pandas series
        abundances['H'] = aux['h1']+2*aux['h2']+3*aux['h3']+4*aux['h4']
        return pd.Series(abundances)

kpar = kppParameters(None)
kpar.read_smiles('mcm_'+root+'_mass.txt')
smiles = pd.Series(kpar.smiles)

sm = SmilesMatch()
stoich = OrderedDict()
for i in range(smiles.shape[0]):
    spec = smiles.index[i].strip()
    a = sm.get_abundances(smiles[i].strip())
    nzv = a > 0
    # formula = ' + '.join(map('{:d}{:s}'.format,a[nzv],a.index[nzv]))
    formula = ' + '.join(map('{}{}'.format,a[nzv].astype(str).replace('1',''),a.index[nzv]))
    stoich[spec] = ('{spec} = {formula} ;'.format(spec=spec,formula=formula))

###_* -------------------- edit kpp file --------------------

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def fix_exponent(formula):
    left = formula.find('**-')+3
    for i in range(left,len(formula)):
        if not is_number(formula[i+1]):
            if formula[i+1] == '.':
                pass
            else:
                right = i
                break
    return formula[0:left-1]+'('+formula[left-1:right+1]+')'+formula[right+1:]

## included in inorganic.kpp
inorg = ['N2O5','H2O2','H2','NA','HONO','SO2','O','HNO3','SO3','O1D','HO2NO2',
         'CO','SA','HSO3','O3','NO','OH','HO2','NO2','NO3','H2O','O2','CL']

kppfile = root+'.kpp'
with open(kppfile,'r') as finp:
    with open(kppfile+'~','w') as fout:
        ## --- header ---
        for line in finp:
            line = line.replace('\r\n','\n')
            if '#INCLUDE atoms' in line:
                line = '{#INCLUDE atoms}' + '\n'
                fout.write(line)
                break
            fout.write(line)
        ## --- add species stoichiometry ---            
        for line in finp:
            line = line.replace('\r\n','\n')
            if '#INLINE F90_RCONST' in line:
                fout.write(line)
                break
            species = line.split()
            if species[0] in inorg:
                continue
            if species[0] in stoich.keys():
                line = stoich[species[0]] + '\n'
            # elif species[0] == 'FPPN': # ??
            #     line = 'FPPN = 10**(LOG10(FCPPN)/(1+(LOG10(KRPPN)/NCPPN)**2))'+'\n'
            fout.write(line)
        ## --- fix exponents ---
        for line in finp:
            line = line.replace('\r\n','\n')
            if '**-' in line:
                line = fix_exponent(line)
            fout.write(line)
os.rename(kppfile+'~',kppfile)
