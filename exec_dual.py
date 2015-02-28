#!/usr/bin/env python

################################################################################
## runs executables
##
## ~exec_dual.py~
##
## USAGE:
## $ python /path/to/exec_dual.py {ROOT} {RUNPATH} {(MODE)}
##
## EXAMPLES:
## $ python ~/git/projects/kppaermod/exec_dual.py apinene 1,2
##
## satoshi.takahama@epfl.ch
##
################################################################################

###_* -------------------- import libraries --------------------

import os
import sys
import subprocess
from collections import OrderedDict
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='run kpp for gas-phase and total simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('RUNPATH',type=str)
parser.add_argument('MODE',type=str,default='gas,total')
args = dict(vars(parser.parse_args()))

###_* -------------------- define paths --------------------

args['MAINPATH'] = os.path.dirname(__file__)
HERE = os.getcwd()
paths = OrderedDict([(p, os.path.join(HERE,'exec_'+p))  for p in ['gas','total']])

## modify environment PATH
env = os.environ.copy()
env['PATH'] = '{}:{}'.format(os.path.join(args['MAINPATH'],'postprocess'),env['PATH'])

###_* -------------------- execute and postprocess --------------------

files = {
    'required':[
        'photolysis.txt'
        ],
    'required_total':[
        'input_partitioning.txt',
        'molefrac_init.txt'
        ],
    'optional':[
        'input_time.txt',
        'input_temp.txt',
        'cgas_init.txt'
        ],
    'outputs':[
        'output_CGAS.txt',
        'output_CAER.txt',
        'output_ERRORS.txt'
        ],
    'dat':[
        '{ROOT}.dat'.format(**args),
        '{ROOT}_aer.dat'.format(**args)
        ]
    }

def ln(src,dst):
    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(src,dst)

def constructfiles(srcpath,dstpath,f):
    return os.path.join(srcpath,f), os.path.join(dstpath,f)

runpath = os.path.join(HERE,args['RUNPATH'])

for label,p in paths.items():
    ## option
    if label not in args['MODE']:
        continue
    ## setup
    execpath = os.path.join(HERE,p)
    outpath = os.path.join(runpath,label)
    linked = []
    ## --- create symlink to inputs ---
    ## required for all
    for f in files['required']:
        src = os.path.join(runpath,f)
        dst = os.path.join(execpath,f)
        if not os.path.exists(src):
            sys.exit('missing required input file: '+f)
        ln(src,dst)
        linked.append(f)
    ## required for total
    if 'total' in args['MODE']:
        for f in files['required_total']:
            src, dst = constructfiles(runpath,execpath,f)
            if not os.path.exists(src):
                sys.exit('missing required input file: '+f)
            ln(src,dst)
            linked.append(f)
    ## 
    for f in files['optional']:
        src, dst = constructfiles(runpath,execpath,f)        
        if os.path.exists(src):
            ln(src,dst)
            linked.append(f)
    ## change directory
    os.chdir(p)
    ## execute run
    subprocess.call('./{ROOT}.exe'.format(**args), shell=True)
    ## move output
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for f in files['outputs']+files['dat']:
        if not os.path.exists(f):
            continue
        dst = os.path.join(outpath,f)
        os.rename(f,dst)
        if f in files['dat']:
            # create formatted file
            subprocess.call('format_output.py '+dst, shell=True, env=env)
    ## clean 
    for f in linked:
        os.remove(f)
