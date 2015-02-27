#!/usr/bin/env python

################################################################################
##
## search_struct.py
##
## run in directory containing mcm_{ROOT}_mass.txt
##
## satoshi.takahama@epfl.ch
##
################################################################################

###_* -------------------- import libraries --------------------

import os
import sys
scriptspath = os.path.join(os.path.dirname(__file__),'scripts')
sys.path.append(scriptspath)
from kpp_generate_SIMPOLGroups import kppParameters
import subprocess
import argparse
import pandas as pd

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='substructure search and vapor pressure calc')
parser.add_argument('ROOT',type=str)
parser.add_argument('PROGPATH',type=str,default='')
args = dict(vars(parser.parse_args()))

args['SCRIPTSPATH'] = scriptspath
args['TMPFILE'] = 'tmp_SMILES.csv'
args['SEARCHEXE'] = os.path.join(args['PROGPATH'],'substructure_search.py')


###_* -------------------- fragments --------------------

kpar = kppParameters(None)
kpar.read_smiles_table('mcm_{ROOT}_mass.txt'.format(**args))

kpar.smiles['SMILES'].reset_index().to_csv(args['TMPFILE'],index=False)
subprocess.call('python {SEARCHEXE} -d -g SIMPOLgroups.csv -i {TMPFILE} -o {ROOT}_SIMPOLGroups.csv'.format(**args), shell=True)
subprocess.call('python {SEARCHEXE} -d -g FTIRgroups.csv -i {TMPFILE} -o {ROOT}_FTIRGroups.csv'.format(**args), shell=True)
os.remove(args['TMPFILE'])

## ###_* -------------------- merge --------------------
# frags = pd.read_csv('{ROOT}_SIMPOLfrags.csv'.format(**args))
# full = pd.merge(kpar.smiles['SMILES'].reset_index(),frags,on='compound')
# full.to_csv('{ROOT}_SIMPOLfrags.csv'.format(**args),index=False)

###_* -------------------- properties --------------------

props = {
    'props_298.csv':298.15,
    'props_358.csv':298.15+60,
}

for filename,temp in props.items():
    sch  = {'OUTFILE':filename,'TEMP':temp}
    sch.update(args)
    subprocess.call('python {SCRIPTSPATH}/simpol2.py -i {ROOT}_SIMPOLGroups.csv -o {ROOT}_{OUTFILE} -t {TEMP:.2f}'.format(**sch),shell=True)
