#!/usr/bin/env python

################################################################################
##
## volatility_bin_matrix.py
##
## run in directory containing mcm_{ROOT}_mass.txt
##
## satoshi.takahama@epfl.ch
##
################################################################################

###_* -------------------- import libraries --------------------

import numpy as np
import pandas as pd
import argparse
from simpol2 import SimpolClass

###_* -------------------- define functions --------------------

def ppb2massc(ppb, molwt, temp=298.15, press=1):
    R = 0.08206
    return ppb*molwt/R*press/temp

def atm2ppb(p):
    ppb = p*1e9
    return ppb

def atm2massc(atm, molwt, temp=298.15, press=1):
    R = 8.206e-5
    return atm*molwt/R*press/temp

def create_binmatrix(c0,bmin=-5,bmax=11,bwidth=1):
    ##
    logbins = np.arange(bmin,bmax+1,bwidth)
    bwidth = float(bwidth)    
    ##
    binmatrix = np.zeros((len(c0),len(logbins)),dtype=int)
    for i in range(len(logbins)):
        ismemberp = np.logical_and(c0 >  10**(logbins[i]-bwidth/2),
                                   c0 <= 10**(logbins[i]+bwidth/2))
        binmatrix[:,i] = np.where(ismemberp,1,0)
    #columns = map('log(C0)={}'.format,logbins))
    table = pd.DataFrame(binmatrix,index=c0.index,columns=logbins)
    return table

if __name__=='__main__':

###_* -------------------- accept arguments --------------------

    parser = argparse.ArgumentParser(description='''
Volatility bin calculation. Provided with abundances of SIMPOL.1 functional groups and molecular weight of compounds, produces an output file denoting which volatility bin in the range log(C^0)=-5,-4,...,11. Example:

$ python ~/git/projects/kppaermod/util/volatility_bin_matrix.py \
compounds/apinene_Sax/mcm_apinene_mass.txt \
compounds/apinene_Sax/apinene_SIMPOLGroups.csv \
outfile.csv

''')
    parser.add_argument('MOLWTFILE',type=str,help='file containing molecular weights, currently reads in MCM format (usually goes by the name "mcm_{root}_mass.txt")')
    parser.add_argument('FRAGTABLE',type=str,help='SIMPOL.1 fragment table')
    parser.add_argument('OUTFILE',type=str,help='name of output file (CSV format)')
    parser.add_argument('TEMP',type=str,nargs='?',default='298.15',help='temperature in Kelvin')
    args = dict(vars(parser.parse_args()))

###_* -------------------- read/calc/export --------------------

    ## read
    simp = SimpolClass()
    simp.read_compounds(args['FRAGTABLE'])
    vp = simp.calc_properties(float(args['TEMP']))['p0']
    ##
    molwt = pd.read_table(args['MOLWTFILE'],
                          skiprows=18,header=None,skipinitialspace=True,
                          names=['compound','SMILES','InChI','molwt'],
                          index_col='compound').molwt
    molwt.index = molwt.index.map(str.strip)

    ## calc
    c0 = ppb2massc(atm2ppb(vp),molwt,float(args['TEMP']))
    binmatrix = create_binmatrix(c0)

    ## export
    binmatrix.reset_index().to_csv(args['OUTFILE'],index=False)
    print '...', args['OUTFILE'], 'created'
