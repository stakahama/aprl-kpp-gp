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

def create_binmatrix(c0,bmin=-2,bmax=11,bwidth=1):
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

    parser = argparse.ArgumentParser(description='(volatility) bin matrix calculation')
    parser.add_argument('ROOT',type=str)
    parser.add_argument('FRAGTABLE',type=str)
    parser.add_argument('OUTFILE',type=str)
    parser.add_argument('TEMP',type=str,nargs='?',default='298.15')
    args = dict(vars(parser.parse_args()))

###_* -------------------- read/calc/export --------------------

    ## read
    simp = SimpolClass()
    simp.read_compounds(args['FRAGTABLE'])
    vp = simp.calc_properties(float(args['TEMP']))['p0']
    ##
    molwt = pd.read_table('mcm_{ROOT}_mass.txt'.format(**args),
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
