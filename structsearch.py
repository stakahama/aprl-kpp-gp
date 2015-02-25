
import os
import sys
sys.path.append(os.path.dirname(__file__),'scripts')
from kpp_generate_SIMPOLGroups.py import kppParameters
import subprocess
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='build kpp for gas-phase and total simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('SEARCH',type=str)
args = dict(vars(parser.parse_args()))

kpar = kppParameters
kpar.read_smiles_table('mcm_{ROOT}_mass.txt'.format(**args))

tmpfile = 'tmp_SMILES.csv'
kpar.smiles['SMILES'].to_csv(tmpfile)
subprocess.call('python {SEARCH} -d -g SIMPOLgroups.csv -i tmp_SMILES.csv -o fragments.csv'.format(**args), shell=True)
os.remove(tmpfile)
