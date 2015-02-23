# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 16:49:52 2014

@author: shipley
"""

import pandas as pd
import numpy as np
import pybel
from collections import OrderedDict

root = raw_input("ROOT: ")


table = pd.read_table('mcm_'+root+'_mass.txt',skiprows=17,header=0,
                      names=['species','SMILES','InChI','molwt'])
                      
smi = []
spec = []
for i in range(0,len(table.species)):
    item1 = table.SMILES[i].strip(' ')
    smi.append(item1)
    item2 = table.species[i].strip(' ')
    spec.append(item2)
    

class Simpolclass:
    def __init__(self):
        ## SMARTS patterns required by SIMPOL
        self.smartspatt = OrderedDict([
            ('C','[#6]'),
            ('O','[#8]'),
            ('H','nomatch'),
            ('N','[#7]'),
            ('S','[#16]'),
            ])
            
        self.smartsaux = OrderedDict([
            ('h1','[*h1]'),
            ('h2','[*h2]'),
            ('h3','[*h3]'),
            ('h4','[*h4]'),
            ])            

    def get_groupnames(self):
	return self.table['groups'].tolist()

    def get_groups(self,smilesstr=None):
        if not smilesstr:
            return np.array([(1 if '_zero'==p else 0) for p in self.smartspatt.keys()])
        ## main body
        mol = pybel.readstring('smi',smilesstr)
        ## store SIMPOL patterns in ordered dictionary
        abundances = OrderedDict()
        for key,patt in self.smartspatt.items()[0:]:
            abundances[key] = 0 if 'nomatch' in patt else \
                              len(pybel.Smarts(patt).findall(mol))
        ## find auxiliary patterns
        aux = OrderedDict()
        for key,patt in self.smartsaux.items():
            aux[key] = len(pybel.Smarts(patt).findall(mol))                      
                              
        abundances['H'] = aux['h1']+2*aux['h2']+3*aux['h3']+4*aux['h4'] 
                              
        return np.array(abundances.values())
    
code = ['C','O','H','N','S']
kpp = []
simp = Simpolclass()
for i in range(0,len(smi)):
    formula = ''
    results = simp.get_groups(smi[i])
    if results[0] == 0:
        pass
    elif results[0] == 1:
        formula += code[0]
    elif results[0] > 1:
        formula += str(results[0]) + code[0]    
    for j in range(1,5):
        if results[j] == 0:
            pass
        elif results[j] == 1:
            formula += ' + '
            formula += code[j]
        else:
            formula += ' + '
            formula += str(results[j]) + code[j]
    
    #print formula
    
    kpp.append(spec[i] +  ' = ' + formula + ' ;')

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

inorg = ['N2O5','H2O2','H2','NA','HONO','SO2','O','HNO3','SO3','O1D','HO2NO2',
         'CO','SA','HSO3','O3','NO','OH','HO2','NO2','NO3','H2O','O2','CL']
    
with open(root+'.kpp','r') as x:
    lines = x.readlines()
    
with open(root+'.kpp','w') as x:
    for line in lines:
        species = line.split()
        if species[0] not in inorg:
            for ln in kpp:
                specieskpp = ln.split()
                if species[0] == specieskpp[0]:
                    line = ln + '\n'
            if species[0] == 'FPPN':
                line = 'FPPN = 10**(LOG10(FCPPN)/(1+(LOG10(KRPPN)/NCPPN)**2))'+'\n'
            elif '#INCLUDE atoms' in line:
                line = '{#INCLUDE atoms}' + '\n'
            elif '**-' in line:
                left = line.find('**-')+3 
                for i in range(left,len(line)):
                    if not is_number(line[i+1]):
                        if line[i+1] == '.':
                            pass
                        else:
                            right = i
                            break
                line = line[0:left-1]+'('+line[left-1:right+1]+')'+line[right+1:]              
            x.write(line)