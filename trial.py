# import argparse

# parser = argparse.ArgumentParser(description='build kpp for gas-phase simulation')
# parser.add_argument('ROOT',type=str)
# parser.add_argument('KPPC',type=str)
# parser.add_argument('CPATH',type=str)
# parser.add_argument('--skipbuild', dest='skipbuild', action='store_true')
# parser.set_defaults(skipbuild=False)
# args = dict(vars(parser.parse_args()))

# print args

import os
import sys
import subprocess

args = {}
args['MAINPATH'] = os.path.dirname(__file__)
scriptspath = os.path.join(args['MAINPATH'],'scripts')
if scriptspath not in sys.path:
    sys.path.insert(0,scriptspath)

env = os.environ.copy()
env['PATH'] = '{}:{}'.format(scriptspath,env['PATH'])

subprocess.call('kpp_extract_indices.py --help',shell=True,env=env)
