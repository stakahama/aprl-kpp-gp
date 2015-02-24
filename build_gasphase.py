#!/usr/bin/env python

import subprocess
import os
import glob
import shutil
import argparse

parser = argparse.ArgumentParser(description='build kpp for gas-phase simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('KPPC',type=str)
parser.add_argument('CPATH',type=str)
args = dict(vars(parser.parse_args()))

args['MAINPATH'] = os.path.dirname(__file__)

# globals().update(args)
# exit()

HERE = os.getcwd()
paths = {p: os.path.join(HERE,p)  for p in ['gas','total']}

for p in paths:
    shutil.rmtree(p)

## runs
runpaths = glob.glob('run_*')

## compounds
cpath = '{CPATH}/{ROOT}'.format(**args)
compoundfiles = [f for f in os.listdir(cpath) if '.'!=f[0]]
subprocess.call('cp -pv {}/* .'.format(cpath), shell=True)

## KPP
os.mkdir(paths['gas'])
for f in compoundfiles:
    os.symlink(os.path.join(HERE,f),os.path.join(paths['gas'],f))
for f in runpaths:
    os.symlink(os.path.join(HERE,f),os.path.join(paths['gas'],f))
    
os.chdir(paths['gas'])
# subprocess.call('cp -pv {MAINPATH}/modules_generic/* .'.format(**args), shell=True)
subprocess.call('cp -pv {MAINPATH}/modules/* .'.format(**args), shell=True)
subprocess.call('cp -pv ../{KPPC} kpp_constants.f90'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/editkpp_v3.py {ROOT}'.format(**args), shell=True)
# subprocess.call('kpp {ROOT}.def'.format(**args), shell=True)

shutil.copytree(paths['gas'],paths['total'],True)

exit()
###_* -------------------- gas only --------------------
subprocess.call('python {MAINPATH}/scripts/kpp_edit_initialize_gasphase.py {ROOT} '.format(**args), shell=True)
subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True)

###_* -------------------- total --------------------
os.chdir(paths['total'])
subprocess.call('python {MAINPATH}/scripts/kpp_edit_initialize_total.py {ROOT} '.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/extract_org_indices_and_molmasses.py {ROOT} mcm_{ROOT}_mass.txt'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_generate_AERmodules.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_main.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_makefile.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_generate_SIMPOLGroups_extended.py {ROOT} mcm_{ROOT}_mass.txt'.format(**args), shell=True)
subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True)
os.chdir(HERE)

## STEP 3 (I doubt you need all of these files):
# 'cp -pv $MAINPATH/extra/* .'
