#!/usr/bin/env python

################################################################################
##
## exec_dual.py
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

outputs = {
    'gas':'{ROOT}.dat'.format(**args),
    'aer':'{ROOT}_aer.dat'.format(**args)
    }

###_* -------------------- parse numbers --------------------

if args['NUMBERS']=='all':
    pass
else:
    numbers = map(int,args['NUMBERS'].split(','))

###_* -------------------- execute and postprocess --------------------

for i in numbers:
    runpath = 'run_{:03d}'.format(i)
    for label,p in paths.items():
        os.chdir(p)
        print os.getcwd()        
        subprocess.call('./{root}.exe {number}'.format(root=args['ROOT'],number=i), shell=True)
        if not os.path.exists(os.path.join(runpath,label)):
            os.mkdir(os.path.join(runpath,label))
        for f in outputs.values():
            if not os.path.exists(f):
                continue
            dst = os.path.join(runpath,label,f)
            os.rename(f,dst)            
            subprocess.call('python {path}/extra/format_output.py {filename}'.format(path=args['MAINPATH'],filename=dst), shell=True)
