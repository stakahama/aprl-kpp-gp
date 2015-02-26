#!/usr/bin/env python

################################################################################
## builds two versions of kpp:
## 1) gas-phase simulations
## 2) total (gas+aerosol) simulations
##
## ~build_dual.py~
##
## USAGE:
## $ python /path/to/build_dual.py {ROOT} {LIGHT} {CPATH}
## {ROOT} = KPP root
## {LIGHT} = photolysis file {original|dark|light}
## {CPATH} = compound path
##
## EXAMPLE:
## $ python ~/git/projects/kppaermod/build_dual.py apinene \
##          dark ../compounds --skipbuild
##
## satoshi.takahama@epfl.ch
##
################################################################################

###_* -------------------- import libraries --------------------

import os
import shutil
import glob
import subprocess
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='build kpp for gas-phase and total simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('LIGHT',type=str)
parser.add_argument('CPATH',type=str)
parser.add_argument('--skipbuild', dest='skipbuild', action='store_true')
parser.set_defaults(skipbuild=False)
args = dict(vars(parser.parse_args()))

# globals().update(args)
# exit()

###_* -------------------- define and create paths --------------------

args['MAINPATH'] = os.path.dirname(__file__)
args['PHOTOPATH'] = os.path.join(args['MAINPATH'],'photolysisfiles',args['LIGHT'])
HERE = os.getcwd()
paths = {p: os.path.join(HERE,'exec_'+p)  for p in ['gas','total']}
template = os.path.join(HERE,'kppbuild')

for p in paths.values():
    if os.path.exists(p):
        shutil.rmtree(p)

if(not args['skipbuild']):

    ## --- kpp build ---
    if os.path.exists(template):
        shutil.rmtree(template)

###_* -------------------- create symlinks for common inputs -------------------

    ## runs
    runpaths = glob.glob('run_*')

    ## compounds
    cpath = '{CPATH}/{ROOT}/*'.format(**args)
    compoundfiles = map(os.path.basename,glob.glob(cpath))
    subprocess.call('cp -pv {} .'.format(cpath), shell=True)

    os.mkdir(template)
    for f in compoundfiles:
        os.symlink(os.path.join(HERE,f),os.path.join(template,f))
    for f in runpaths:
        os.symlink(os.path.join(HERE,f),os.path.join(template,f))

###_* -------------------- build kpp (common) --------------------

# print '-------------------- building kpp common --------------------'    
    os.chdir(template)
    subprocess.call('cp -pv {MAINPATH}/modules_generic/* .'.format(**args), shell=True)
    subprocess.call('cp -pv {PHOTOPATH}/photolysis.txt .'.format(**args), shell=True)
    subprocess.call('python {MAINPATH}/scripts/kpp_edit_kpp.py {ROOT}'.format(**args), shell=True)
    subprocess.call('kpp {ROOT}.def'.format(**args), shell=True)
    ## ---

shutil.copytree(template,paths['gas'],True)
shutil.copytree(template,paths['total'],True)

###_* -------------------- gasphase only --------------------

# print '-------------------- building gas-phase --------------------'
os.chdir(paths['gas'])
subprocess.call('python {MAINPATH}/scripts/kpp_edit_initialize_gasphase.py {ROOT} '.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_makefile_gasphase.py {ROOT}'.format(**args), shell=True)
subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True)

###_* -------------------- total --------------------

# print '-------------------- building total --------------------'
os.chdir(paths['total'])
subprocess.call('cp -pv {MAINPATH}/modules_partitioning/* .'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_initialize.py {ROOT} '.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_makefile.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_extract_org_indices_and_molmasses.py {ROOT} mcm_{ROOT}_mass.txt'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_generate_AERmodules.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_edit_main.py {ROOT}'.format(**args), shell=True)
subprocess.call('python {MAINPATH}/scripts/kpp_generate_SIMPOLGroups.py {ROOT} mcm_{ROOT}_mass.txt {ROOT}_SIMPOLGroups.csv'.format(**args), shell=True)
subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True)
os.chdir(HERE)

## STEP 3 (I doubt you need all of these files):
# 'cp -pv $MAINPATH/extra/* .'
