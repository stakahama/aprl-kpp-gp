#!/usr/bin/env python

################################################################################
##
## extracts indices of the concentration vectors
## (+molecular weights for organic compounds)
## formats output for csv file and/or f90 file
##
## S. Takahama and F.Bernhard 02.01.15 in order to generate
## additionally molecular masses
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################

###_* -------------------- import libraries --------------------

import os
import sys
from collections import OrderedDict
from kpp_generate_SIMPOLGroups import kppParameters
import pandas as pd
import argparse

###_* -------------------- parse arguments --------------------

parser = argparse.ArgumentParser(description='extracts indices and molecular masses of the concentration vectors')
parser.add_argument('root',type=str,help="used to read in {root}_Parameters.f90")
parser.add_argument('SMILESfile',type=str,nargs='?',default='',help="[optional] name of mcm_subset_mass.txt")
parser.add_argument('mode',type=int,nargs='?',default=0,help="[optional] one of 0: both compound_indices_table.csv and org_indices.f90 (default); 1: compound_indices_table.csv only; 2: org_indices.f90 only")
args = parser.parse_args()
root = args.root
smilesfile = args.SMILESfile
mode = args.mode

###_* -------------------- read {root}_Parameters.f90 --------------------
parms = kppParameters(root)
parms.read_parms()                                  #-> parms.ind
compounds = OrderedDict(map(reversed,parms.ind))    # should be sorted
compounds_table = pd.DataFrame(compounds.items(),columns=['compound','index']).set_index('compound')

###_* -------------------- read and process {SMILESfile} --------------------

if smilesfile != '':
    parms.read_smiles_table(smilesfile)             #-> parms.smiles
    mol_masses = parms.smiles.molwt.values          # molecular masses of organics
    indices = map(compounds.get,parms.smiles.index) # compound indices
    ## add to compounds_table
    compounds_table['organic'] = 0
    compounds_table.ix[parms.smiles.index,'organic'] = 1

###_* -------------------- export --------------------

def write_vector_f90(fout,label,vector,indent,nmaxline):
    nvec = len(vector)
    fout.write(' allocate({:s}({:d}))\n'.format(label,len(vector)))
    fout.write(' {:s}= (/'.format(label))
    st,en = 0,nmaxline
    while st <= nvec:
        if st > 0:
            fout.write(indent)
        fout.write(','.join(map(str,vector[st:en])))
        if en < nvec:
            fout.write(',&\n')
        st,en = en,en+nmaxline
    fout.write('/)'+'\n')

###_ . compounds table

## output of compound names and indices (useful for preparing initial concentrations, mole fractions)
if mode < 2:
    compounds_table.reset_index().to_csv('compound_indices_table.csv',index=False)

###_ . f90 file
if mode in [0,2] and smilesfile != '':
    indent = ' '*7
    with open('org_indices.f90','w') as fout:
        write_vector_f90(fout,'organics',indices,indent,15)
        write_vector_f90(fout,'organics_mol_masses',mol_masses,indent,5)
