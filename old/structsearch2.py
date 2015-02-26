

###_* -------------------- import library --------------------

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),'scripts'))
sys.path.append(os.path.join(os.path.dirname(__file__),'lib'))
from kpp_generate_SIMPOLGroups import kppParameters
import subprocess
import argparse
from simpolmoddebug import Simpolclass

import simpolmoddebug as ss
reload(ss)
###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='build kpp for gas-phase and total simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('SEARCH',type=str)
args = dict(vars(parser.parse_args()))

###_* -------------------- do stuff --------------------

kpar = kppParameters(None)
kpar.read_smiles_table('mcm_{ROOT}_mass.txt'.format(**args))

# tmpfile = 'tmp_SMILES.csv'
# kpar.smiles['SMILES'].reset_index().to_csv(tmpfile,index=False)
# subprocess.call('python {SEARCH} -d -g SIMPOLgroups.csv -i tmp_SMILES.csv -o fragments.csv'.format(**args), shell=True)
# os.remove(tmpfile)

import pandas as pd

simp = ss.Simpolclass()
header = map(lambda x: x[1:],simp.get_groupnames())
frag = kpar.smiles['SMILES'].apply(lambda x: pd.Series(simp.get_groups(x),index=header))
frag.reset_index().to_csv('fragments2.csv',index=False)

