#!/usr/bin/env python

################################################################################
##
## kpp_generate_SIMPOLGroups.py
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
################################################################################

import os
# import sys
# dd = lambda x: os.path.dirname(os.path.dirname(x))
# sys.path.append(os.path.join(dd(__file__),'lib'))
import re
from operator import itemgetter
from collections import OrderedDict
import numpy as np
import pandas as pd
import pybel
import argparse

parser = argparse.ArgumentParser(description='create a {ROOT}_SIMPOLGroups.f90 containing SIMPOL group abundances for each molecule')
parser.add_argument('root',type=str)
parser.add_argument('SMILESfile',type=str)
parser.add_argument('SMARTSfile',type=str)

class kppParameters:
    
    def __init__(self,root,path='.'):
        ## prefix is the value of ROOT used as example in KPP
        ##   manuscript by Sandu
        ## path is for both ROOT_Parameters.f90 file (input)
        ##   and ROOT_SIMPOLGroups.f90
        self.root = root
        self.path = path
        self.outvar = 'out' # used by Gencase and Writemodule
        
    def read_parms(self):
        ## reads {root}_Parameters.f90 file
        ## stores list containing pairs of indices and MCM molecule names
        def extract(line):
            return patt.sub('',line).strip().split(' = ')[::-1]
        def asTuple(pair):
            return int(pair[0]),pair[1].replace('ind_','')
        filename = os.path.join(self.path,self.root+'_Parameters.f90')
        patt = re.compile('.+:: ')
        start = '! Index declaration for variable species in C and VAR'
        with open(filename) as f:
            for line in f:     # read until match
                if start in line:
                    break
            for i in range(2): # skip 2 lines
                next(f)
            ind = []
            for line in f:     # read parameters
                if '::' not in line:
                    break
                ind.append(asTuple(extract(line)))
        self.ind = ind

    def read_smiles(self,smilesfile,attr=('smiles',range(2))):
        ## to change to molwt: set attr=('molwt',[0,3])
        ## smilesfile ('mcm_subset_mass.txt')
        with open(smilesfile) as f:
            ## skip header
            count = 0
            for line in f:
                if '*****' in line:
                    count += 1
                    if count == 3:
                        break
            ## save MCM names and smiles string (first two fields)
            pairs = []
            for line in f:
                if len(line) > 2:
                    pairs.append(map(str.strip,itemgetter(*attr[1])(line.strip().split('\t'))))
        setattr(self,attr[0],OrderedDict(pairs))

    def read_smiles_table(self,smilesfile,attr='smiles',columns=['compound','SMILES','InChI','molwt']):
        import pandas as pd        
        ## smilesfile ('mcm_subset_mass.txt')
        with open(smilesfile) as f:
            ## skip header
            count = 0
            for line in f:
                if '*****' in line:
                    count += 1
                    if count == 3:
                        break
            ## save MCM names and smiles string (first two fields)
            fields = []
            for line in f:
                if len(line) > 2:
                    fields.append(map(str.strip,line.strip().split('\t'))) # convert to tuple(?)
        setattr(self,attr,pd.DataFrame(fields,columns=columns).set_index('compound'))

    def gen_case(self,groupfn):
        ## generates fortran 90 case statement
        ## requires function which converts MCM name to list containing
        ## integer of group abundance
        ## stores case statement as multi-line text string
        def format_out(data,n=16):
            # formats list {data} of 31 integers (moles of SIMPOL groups)
            # returns string as '{outvar}=(/.../)'
            indent = ' '*10
            strlist = map(str,data)
            interv = range(n,len(data),n)
            if len(data) % n > 0:
                interv.append(len(data))
            strout = ''
            start = 0
            for count,end in enumerate(interv):
                if count > 0:
                    strout += indent
                strout += ','.join(strlist[start:end])
                if count < len(interv)-1:
                    strout += ',&\n'
                start = end
            return indent+'{}=(/{}/)'.format(self.outvar,strout)
        ## print case statements for each MCM molecule
        newline = '\n'
        indent = ' '*7
        case = ''
        for p in self.ind:
            smi = self.smiles[p[1]] if p[1] in self.smiles.keys() else None
            case += indent+'case ({:d}) ! {:s}, {:s}'.format(p[0],p[1],smi)+newline
            # case += format_out(groupfn(smi))+newline
            case += format_out(groupfn(p[1]))+newline            
        case += indent+'case default'+newline
        case += format_out([0]*31)
        self.select_case = case
        
    def write_module(self,path=None):
        ## writes a f90 module file with _SIMPOLGroups.f90 extension
        if not path: # debug
            path = self.path
        filename = os.path.join(path,self.root+'_SIMPOLGroups.f90')
        module_template = '''module {ROOT}_SIMPOLGroups
        
  integer, parameter :: ngroups = 31

contains

  function substruct_count(index) result({OUT})
  
    integer, intent(in) :: index
    integer, dimension(ngroups) :: {OUT}
    
    assignvector: select case (index)
{CASE_STATEMENT}
    end select assignvector
    
  end function substruct_count

end module {ROOT}_SIMPOLGroups
'''
        with open(filename,'w') as f:
            f.write(module_template.format(ROOT=self.root,
                                           OUT= self.outvar,
                                           CASE_STATEMENT=self.select_case))
        print filename+' created'

class SMARTSgetter:

    def __init__(self,filename):
        self.table = pd.read_csv(filename,index_col='compound')

    def re_index(self,smilesdict):
        self.table.index = map(smilesdict.get,self.table.index)

    def get_groups(self,compound):
        if compound in self.table.index:
            return self.table.ix[compound].values
        else:
            return np.array([(1 if p=='zero' else 0) for p in self.table.columns])

if __name__=='__main__':

    args = parser.parse_args()
    root = args.root
    ## root, smilesfile = ('octane_gen','mcm_subset_mass.txt')
    parms = kppParameters(root)    
    parms.read_parms()
    parms.read_smiles(args.SMILESfile)
    # parms.smiles = OrderedDict(smarts.table['SMILES'])
    # smarts.table.drop('SMILES',axis=1,inplace=True)
    smarts = SMARTSgetter(args.SMARTSfile)
    # smarts.re_index(parms.smiles)
    parms.gen_case(smarts.get_groups)
    parms.write_module(path='.')
