################################################################################
##
## kpp_generate_AERmodules.py
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
################################################################################

import argparse
import os

parser = argparse.ArgumentParser(description='create several modules for handling additional aerosol variables')
parser.add_argument('root',type=str)

template = lambda x: os.path.join(os.path.dirname(__file__),'templates',x)

class AERModules:

    def __init__(self,root):
        
        self.root = root
        self.modules = [
            'GlobalAER',
            'InitializeAER',
            'UtilAER',
            'PartitionAER'
            ]
    
    def write(self):
        for mod in self.modules:
            with open(template('ROOT_{MODULE}.f90'.format(MODULE=mod))) as finp:
                contents = finp.read().format(ROOT=self.root)
            outfile = '{ROOT}_{MODULE}.f90'.format(ROOT=self.root, MODULE=mod)
            print outfile+' created'
            with open(outfile,'w') as fout:
                fout.write(contents)

if __name__=='__main__':

    args = parser.parse_args()
    root = args.root
    ## root = 'octane_gen'
    mods = AERModules(root)
    mods.write()
