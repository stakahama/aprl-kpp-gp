################################################################################
##
## extracts indices of the concentration vectors
## which correspond to organic compounds only
##
## extended by F.Bernhard 02.01.15 in order to generate
## additionally molecular masses
##
################################################################################

from collections import OrderedDict                                                      # FB: for OrderedDict
#from kpp_generate_SIMPOLGroups import kppParameters
from kpp_generate_SIMPOLGroups_extended import kppParameters
from simpol import Simpolclass
import pandas as pd
import argparse

## root, smilesfile = ('octane_gen','mcm_subset_mass.txt')
#root, smilesfile = ('tm135b','mcm_tm135b_mass.txt')                   # FB: TODO: generate it automatically, similar to "kpp_generate_SIMPOLGroups.py"
parser = argparse.ArgumentParser(description='extracts indices and molecular masses of the concentration vectors')
parser.add_argument('root',type=str)
parser.add_argument('SMILESfile',type=str)
args = parser.parse_args()
root = args.root
smilesfile = args.SMILESfile


## execute
simp = Simpolclass()
# kppParameters.root = root
parms = kppParameters(root)    
parms.read_parms() #-> parms.ind
parms.read_smiles_table(smilesfile) #-> parms.smiles
molar_masses = parms.smiles.molwt.values # FB: generate list with molecular masses of organics

## get compound indices
compounds = OrderedDict(map(reversed,parms.ind))            # FB: needed to import collections
indices = map(compounds.get,parms.smiles.index)

## format output for fortran
indent = ' '*7
with open('org_indices.f90','w') as fout:
    ## indices
    n = 15    
    fout.write(' allocate(organics({:d}))\n'.format(len(indices)))
    fout.write(' organics= (/')
    st,en = 0,n
    while st <= len(indices):
        if st > 0:
            fout.write(indent)
        fout.write(','.join(map(str,indices[st:en])))
        if en < len(indices):
            fout.write(',&\n')
        st,en = en,en+n
    fout.write('/)'+'\n')
    # FB:
    ## format molecular masses output for fortran
    n = 5
    fout.write(' allocate(organics_mol_masses({:d}))\n'.format(len(molar_masses)))        
    fout.write(' organics_mol_masses= (/')
    st,en = 0,n
    while st <= len(molar_masses):
        if st > 0:
            fout.write(indent)
        fout.write(','.join(map(str,molar_masses[st:en])))
        if en < len(molar_masses):
            fout.write(',&\n')
        st,en = en,en+n
    fout.write('/)'+'\n')

## output for preparing initial concentrations (mole fractions)
compounds_table = pd.DataFrame(compounds.items(),columns=['compound','index']).set_index('compound')
compounds_table['organic'] = 0
compounds_table.ix[parms.smiles.index,'organic'] = 1
compounds_table.reset_index().to_csv('org_indices_table.csv',index=False)
