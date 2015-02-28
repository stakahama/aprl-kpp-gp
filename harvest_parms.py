#!/usr/bin/env python

###_* -------------------- import libraries --------------------
import os
import re
import glob
import pandas as pd
from collections import OrderedDict
import argparse

###_* -------------------- parse arguments --------------------
parser = argparse.ArgumentParser(description='generate a0.txt')
parser.add_argument('ROOT',type=str)
parser.add_argument('RUNPATH',type=str,nargs='*')
args = parser.parse_args()

deffile = args.ROOT+'.def'

runs = args.RUNPATH
if len(runs)==0:
    runs = glob.glob('run_*')

###_* -------------------- define paths --------------------

HERE = os.getcwd()
mech = os.path.basename(HERE)

###_* -------------------- input file definitions --------------------
input_files = {
    'input_time.txt':('TSTART','DURATION','DT'),
    'input_temp.txt':('TEMP','CFACTOR'),
    'input_partitioning.txt':('M0','PARTITION_SUBSTEPS')
    }

###_* -------------------- define functions --------------------

# def update(env,x):
#     variable = x[0].strip()
#     try:
#         value = float(x[1])
#     except:
#         value = eval(x[1],None,env)
#     return OrderedDict(env.items()+[(variable,value)])

def expr2pair(x):
    fields = map(str.strip,x.strip().strip(';').split('='))
    return [(fields[0].upper(),fields[1])]

def modify_time_end(parms):
    t_end = parms['TEND']
    pretty = '{:.0f}'.format
    if 'DURATION' in parms.keys():
        t_end = pretty(eval('{TSTART}+{DURATION}'.format(**parms)))
    else:
        try:
            float(t_end)
        except:
            tokens = re.compile('('+'|'.join(parms.keys())+')')
            t_end = pretty(eval(tokens.sub('{\\1}',t_end).format(**parms)))
    return t_end

###_* -------------------- read def file --------------------

parms = OrderedDict()

with open(deffile) as f:
    for line in f:
        if '#INLINE F90_INIT' in line:
            break
    for line in f:
        if '#ENDINLINE' in line:
            break
        # parms = update(parms,line.strip().strip(';').split('='))
        parms.update(expr2pair(line))
    for line in f:
        if '#INITVALUES' in line:
            break
    for line in f:
        # parms = update(parms,line.strip().strip(';').split('='))
        parms.update(expr2pair(line))

###_* -------------------- harvest runs --------------------

master = None

for runpath in runs:

    augmented = OrderedDict()

    with open(os.path.join(runpath,'photolysis.txt')) as f:
        lightmode = f.readline().strip().replace('#','')
    augmented.update({'LIGHTMODE':lightmode})

    fname = os.path.join(runpath,'molefrac_init.txt')
    if os.path.exists(fname):
        with open(fname) as f:
            molefrac_label = f.readline().strip().replace('#','')
        augmented.update({'MOLEFRAC_INIT':molefrac_label})

    fname = os.path.join(runpath,'cgas_init.def')
    if os.path.exists(fname):
        pairs = []
        with open(fname) as f:
            for line in f:
                pairs += expr2pair(line)
        augmented.update(pairs)

    for inp in input_files.keys():
        fname = os.path.join(runpath,inp)
        if os.path.exists(fname):
            with open(fname) as f:
                pairs = OrderedDict(zip(input_files[inp],f.read().split()))
            augmented.update(pairs)

    parmscopy = parms.copy()
    for k in augmented.keys():
        parmscopy[k] = augmented[k]

    parmscopy['TEND'] = modify_time_end(parmscopy)
    parmscopy['DURATION'] = None

    parmstable = pd.DataFrame(parmscopy.items(),columns=['PARAMETER','VALUE'])
    parmstable.insert(0,'RUN',runpath)
    parmstable.insert(0,'MECHANISM',mech)

    if master:
        master = master.join(parmstable,how='outer')
    else:
        master = parmstable

###_* -------------------- export --------------------

master.to_csv('parameter_table.csv',index=False)
