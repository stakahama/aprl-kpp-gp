#!/usr/bin/env python

################################################################################
## runs executables
##
## ~exec_dual.py~
##
## USAGE:
## $ python /path/to/exec_dual.py {ROOT} {NUMBERS}
## {NUMBERS} can be a single number or several concatenated by commas.
##
## EXAMPLES:
## $ python ~/git/projects/kppaermod/exec_dual.py apinene 1,2
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
        ## create symlink to inputs
        if not os.path.exists(os.path.join(HERE,p,runpath)):
            os.symlink(os.path.join(HERE,runpath),
                       os.path.join(HERE,p,runpath))
        ## change directory
        os.chdir(p)
        ## execute run
        subprocess.call('./{root}.exe {number}'.format(root=args['ROOT'],number=i), shell=True)
        ## move output and created formatted file
        if not os.path.exists(os.path.join(runpath,label)):
            os.mkdir(os.path.join(runpath,label))
        for f in outputs.values():
            if not os.path.exists(f):
                continue
            dst = os.path.join(runpath,label,f)
            os.rename(f,dst)            
            subprocess.call('python {path}/postprocess/format_output.py {filename}'.format(path=args['MAINPATH'],filename=dst), shell=True)
