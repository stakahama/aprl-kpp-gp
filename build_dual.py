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
import sys
import shutil
import glob
import subprocess
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='build kpp for gas-phase and total simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('CPATH',type=str)
parser.add_argument('--skipbuild', dest='skipbuild', action='store_true')
parser.add_argument('--onlygas', dest='onlygas', action='store_true')
parser.add_argument('--onlytotal', dest='onlytotal', action='store_true')
parser.set_defaults(skipbuild=False)
parser.set_defaults(onlygas=False)
parser.set_defaults(onlytotal=False)
args = dict(vars(parser.parse_args()))

# globals().update(args)
# exit()

###_* -------------------- define and create paths --------------------

args['MAINPATH'] = os.path.dirname(__file__)
HERE = os.getcwd()
paths = {p: os.path.join(HERE,'exec_'+p)  for p in ['gas','total']}
template = os.path.join(HERE,'kppbuild')

## modify environment PATH
env = os.environ.copy()
env['PATH'] = '{}:{}'.format(os.path.join(args['MAINPATH'],'scripts'),env['PATH'])

## prepare directories
# for p in paths.values():
#     if os.path.exists(p):
#         shutil.rmtree(p)

if not (args['skipbuild'] or args['onlygas'] or args['onlytotal']):

    ## --- kpp build ---
    if os.path.exists(template):
        shutil.rmtree(template)

###_* -------------------- build def file -------------------

    deffile = '{ROOT}.def'.format(**args)
    
    with open(os.path.join(HERE,deffile),'w') as fout:
        with open(os.path.join(args['MAINPATH'],'scripts','templates','generic.def')) as finp:
            fout.write(finp.read().format(**args))

###_* -------------------- create symlinks for common inputs -------------------

    ## runs
    # runpaths = glob.glob('run_*')

    ## compounds
    compoundfiles = [deffile]
    for f in os.listdir(args['CPATH']):
        if f!=deffile:
            compoundfiles.append(f)
            shutil.copy2(os.path.join(args['CPATH'],f),'.')

    os.mkdir(template)
    for f in compoundfiles:
        os.symlink(os.path.join(HERE,f),os.path.join(template,f))
    # for f in runpaths:
    #     os.symlink(os.path.join(HERE,f),os.path.join(template,f))

###_* -------------------- build kpp (common) --------------------

# print '-------------------- building kpp common --------------------'    
    os.chdir(template)
    subprocess.call('cp -pv {MAINPATH}/modules_generic/* .'.format(**args), shell=True)
    subprocess.call('kpp_edit_kpp.py {ROOT}'.format(**args), shell=True, env=env)
    subprocess.call('kpp {ROOT}.def'.format(**args), shell=True)
    ## ---

    ## test build
    if '{ROOT}_Parameters.f90'.format(**args) not in os.listdir('.'):
        sys.exit('KPP build fail')

    ## indices
    if not os.path.exists(os.path.join(HERE,'compound_indices_table.csv')):
        subprocess.call('kpp_extract_indices.py {ROOT} mcm_{ROOT}_mass.txt 1'.format(**args), shell=True, env=env)
        shutil.move('compound_indices_table.csv',HERE)


###_* -------------------- gasphase only --------------------

# print '-------------------- building gas-phase --------------------'
if not args['onlytotal']:
    if os.path.exists(paths['gas']):
        shutil.rmtree(paths['gas'])
    shutil.copytree(template,paths['gas'],True)
    os.chdir(paths['gas'])
    subprocess.call('kpp_edit_initialize.py {ROOT} gas'.format(**args), shell=True, env=env)
    subprocess.call('kpp_edit_makefile_gasphase.py {ROOT}'.format(**args), shell=True, env=env)
    subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True)

###_* -------------------- total --------------------

# print '-------------------- building total --------------------'
if not args['onlygas']:
    if os.path.exists(paths['total']):
        shutil.rmtree(paths['total'])
    shutil.copytree(template,paths['total'],True)
    os.chdir(paths['total'])
    subprocess.call('cp -pv {MAINPATH}/modules_partitioning/* .'.format(**args), shell=True)
    subprocess.call('kpp_edit_initialize.py {ROOT} total'.format(**args), shell=True, env=env)
    subprocess.call('kpp_edit_makefile.py {ROOT}'.format(**args), shell=True, env=env)
    subprocess.call('kpp_extract_indices.py {ROOT} mcm_{ROOT}_mass.txt 2'.format(**args), shell=True, env=env)
    subprocess.call('kpp_generate_AERmodules.py {ROOT}'.format(**args), shell=True, env=env)
    subprocess.call('kpp_edit_main.py {ROOT}'.format(**args), shell=True, env=env)
    subprocess.call('kpp_generate_SIMPOLGroups.py {ROOT} mcm_{ROOT}_mass.txt {ROOT}_SIMPOLGroups.csv'.format(**args), shell=True, env=env)
    subprocess.call('make -f Makefile_{ROOT}'.format(**args), shell=True, env=env)
    os.chdir(HERE)

## STEP 3 (I doubt you need all of these files):
# 'cp -pv $MAINPATH/extra/* .'
