#!/usr/bin/env python

################################################################################
## runs executables
##
## ~reinitialize_a0.py~
##
## 
##
## satoshi.takahama@epfl.ch
##
################################################################################


###_* -------------------- import libraries --------------------

import os
import sys
scriptspath = os.path.join(os.path.dirname(os.path.dirname(__file__)),'scripts')
sys.path.append(scriptspath)
from simpol2 import SimpolClass
from collections import OrderedDict
import subprocess
import pandas as pd
import numpy as np
import argparse

###_* -------------------- accept arguments --------------------

parser = argparse.ArgumentParser(description='generate a0.txt')
parser.add_argument('ROOT',type=str)
parser.add_argument('RUNPATH',type=str)
args = dict(vars(parser.parse_args()))

###_* -------------------- define paths --------------------

a0file = 'molefrac_init.txt'

args['MAINPATH'] = os.path.dirname(__file__)
HERE = os.getcwd()
paths = OrderedDict([(p, os.path.join(HERE,'exec_'+p))  for p in ['gas','total']])

###_* -------------------- calculate partitioning --------------------

exit()
