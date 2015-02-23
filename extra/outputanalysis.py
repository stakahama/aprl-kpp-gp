# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 12:22:31 2014

@author: shipley
"""

import pandas as pd
import numpy as np
from matmul import matmul
import matplotlib.pyplot as plt
import math
import os.path
from simpolmod import Simpolclass
import pylab
simp = Simpolclass()

root = raw_input("ROOT: ")
phase = raw_input('aer or gas phase? ')

if phase == 'aer':
    phaseresp = 'aer_'
elif phase == 'gas':
    phaseresp = ''
else:
    print 'invalid response!'
    root = 'none'
    
    
keyword = root+'_'+phaseresp

table = pd.read_table('mcm_'+root+'_mass.txt',skiprows=17,header=0,
                      names=['species','SMILES','InChI','molwt'])

inorg = ['N2O5','H2O2','H2','NA','HONO','SO2','O','HNO3','SO3','O1D','HO2NO2',
         'CO','SA','HSO3','O3','NO','OH','HO2','NO2','NO3','H2O','O2','CL']

spec = []                 
smi = []
molwt = []
data1 = {}
for i in range(0,len(table.species)):
    item1 = table.species[i].strip(' ')
    item2 = table.SMILES[i].strip(' ')
    item3 = table.molwt[i]
    if item1 not in inorg:
        spec.append(item1)
        smi.append(item2)
        molwt.append(item3)
        #data1[item1]=item2
    else:
        print 'error in r1', item1
        

###############################################



data1 = {'Species': spec, 'SMILES': smi}

results1 = pd.DataFrame(data1,columns=['Species', 'SMILES'])

with open(root+'_SMILES.txt','w') as r1:
    results1.to_csv(r1,index=True)


out = []
with open (keyword+'formatted.txt', 'r') as r3:
    for line in r3:
        MCM = line.strip().split(',')
        MCMnew = [item for item in MCM if item not in inorg]
        break
    MCMnew.pop(0)

    sminew = []
    for i in MCMnew:
        if i != "TIME":
            num = spec.index(i)
            sminew.append(smi[num])
    
datax = pd.DataFrame.from_csv(keyword+'formatted.txt', header=0, index_col=0)

head = list(datax.columns.values)

timeind = list(datax.index.values)

drop =[z for z in head if z in inorg]            
datax.drop(drop, 1, inplace = True)            

newhead = list(datax.columns.values)

for i in range(0,len(newhead)):
    if newhead[i] in MCMnew:
        ind = MCMnew.index(newhead[i])
        newhead[i]=sminew[ind]

datax.columns = newhead


with open(keyword+'edited.txt','w') as s1:
    datax.to_csv(s1,index=True)
    
newmolwt=[] 
for item in newhead:
    indx = smi.index(item)
    newmolwt.append(molwt[indx] / 0.08206 / 298)

wtfactors = pd.DataFrame(newmolwt,index=newhead,columns=['molwt'])

with open(keyword+'molwts.txt','w') as s2:
    wtfactors.to_csv(s2,index=True)

orgmass=[]    
orgmass = matmul(keyword+'edited.txt',keyword+'molwts.txt')

orgmassdf = pd.DataFrame(orgmass,index=timeind,
                         columns=['total mass in '+phaseresp+' phase'])

    
data2 = []
for i in newhead:
    line = [i]
    line.extend(simp.calc_properties(i,298.15))
    line.extend(simp.get_groups(i))
    line.extend(simp.get_FTIRgroups(i))    
    data2.append(line)

header = ['SMILES','p0(atm)','deltaH(kJ/mol)','zero', 'carbon number',
'carbon number on the acid-side of an amide (asa)','aromatic ring',
'non-aromatic ring', 'C=C (non-aromatic)', 'C=C-C=O in non-aromatic ring',
'hydroxyl (alkyl)', 'aldehyde', 'ketone', 'carboxylic_acid', '_ester',
'_ether', 'ether (alicyclic)', 'ether, aromatic', 'nitrate', 'nitro',
'aromatic hydroxyl', 'amine, primary', 'amine, secondary', 'amine, tertiary',
'amine, aromatic', 'amide, primary', 'amide, secondary', 'amide, tertiary',
'carbonylperoxynitrate', 'peroxide', 'hydroperoxide', 'carbonylperoxyacid',
'nitrophenol', 'nitroester','ftirC','ftirO','ftirH','ftirsecondary amine', 
'ftirtertiary amine','ftirCH2','ftirCH3','ftirtertiary C','ftirprimary amine','ftiralkane CH',
'ftiralkene CH','ftiraromatic CH','ftircarbonyl','ftiralcohol','ftircarboxylic acid','ftirester','ftirether']

results2 = pd.DataFrame(data2,columns=header)

FTIR = []
for i in range(len(data2)):
    FTIR.append(data2[i][34:49])


with open(root+'_Groups.txt','w') as r2:
    results2.to_csv(r2,index=None)
    

results3 = matmul(keyword+'edited.txt',root+'_Groups.txt')
#results2 = matmul('test_aer.txt','COCount.txt')



############### GENERATE FINAL RESULTS TABLE AND FILE ######################
#NOTE: Simpol groups are by molar basis (ppb), FTIR groups are by mass basis (ug/m^3)

### eliminate 'SMILES' label from header, make DF from data, replace index label with 'time'
header.pop(0)
total = pd.DataFrame(results3,index=timeind,columns=header) #+ results2
total.index.name = 'Time(hrs)'

### Generate Results for O:C and H:C Plots
OC = total.ftirO / total.ftirC
HC = total.ftirH / total.ftirC
OCHC = pd.DataFrame(OC,columns=['OC'])
OCHC['HC']=HC

massfactors = {'ftirC':12,'ftirO':16,'ftirH':1,'ftiralkene CH':13,'ftiraromatic CH':8.7,
               'ftirprimary amine':11,'ftiralkane CH':7,'ftircarbonyl':28,'ftiralcohol':23,
               'ftircarboxylic acid':45,'ftirester':45,'ftirether':40}

totalmass=pd.DataFrame(index=total.index)
for i in total.columns.values:
    if i in massfactors:
        factor = massfactors[i]/298.15/0.08206
        mass=pd.DataFrame(total[i]*factor,columns=[i])
        #break
        totalmass[i]=mass


### put results into _RESULTS.txt file (csv)
with open(keyword+'FTIRRESULTS.txt','w') as r3:
    total.to_csv(r3,index=True)

############### PLOT STACKED BAR CHART WITH FTIR GROUPS ####################
plt.close('all')

plotdf = totalmass.ix[:,3:]

colors = ['g','b','c','y',[1,0,.7],'r',[0,1,0],[1,.6,.2],
         'm',[.8,.8,.3],[.5,.5,.1]]
p1 = plotdf.plot(title=keyword,legend=True,kind='bar',stacked=True,color=colors)
plt.legend(loc=8,ncol=3)
plt.ylabel('ug/m^3')

#plt.show()
plt.savefig(keyword+'FTIRfigure.pdf')

    
############### PLOT LINE GRAPH WITH O:C AND H:C RATIOS ####################
#plt.close()

OCHC.plot(title='O:C and H:C atomic ratios for '+keyword,
          legend=True,kind='line',color='br')
plt.legend(loc=5)
plt.ylabel('Atomic Ratio')

plt.savefig(keyword+'OCHC_v_time.pdf')

OCHC.plot(x='OC', y='HC',title='O:C and H:C atomic ratios for '+keyword,
          legend=True,kind='scatter',color='b')

#plt.show()

plt.savefig(keyword+'HC_v_OC.pdf')


orgmassdf.plot(x=orgmassdf.index,title='OM in '+root+' '+phaseresp[0:3]+' phase',
          legend=False,kind='line',color='k')
plt.ylabel('OM (ug/m^3)')
plt.xlabel('Time (hrs)')          

plt.show()

plt.savefig(keyword+'OM.pdf')

#############Make saturation concentration distribution data#################

newarray = datax.tail(n=1).transpose()
newarray.columns = ['ppb']
logC = []
conc = []
j=0
for i in newarray.index:
    find = smi.index(i)
    p0 = simp.calc_properties(i,298)[0]
    logC.append(math.log10(molwt[find]*p0*(10**6)/0.00008206/298))
    ppb = newarray.ix[j,'ppb']
    conc.append(molwt[find]*ppb/0.08206/298)
    j+=1
newarray['log C0'] = logC
newarray['conc (ug/m^3)'] = conc

with open(keyword+'satconcdist.txt','w') as r4:
    newarray.to_csv(r4,index=True)

#############Make saturation concentration distribution plot#################

histlist = [0,0,0,0,0,0,0,0,0,0]
#histdict = {-3:0,'-2':0,'-1':0,'0':0,'1':0,'2':0,'3':0,'4':0,'5':0,'6':0}
for i in range(len(newarray)):
    C = newarray.ix[i,'log C0']
    if C > -20 and C <-2.5:
        histlist[0] += newarray.ix[i,'conc (ug/m^3)']
    elif C <-1.5:
        histlist[1] += newarray.ix[i,'conc (ug/m^3)']
    elif C <-0.5:
        histlist[2] += newarray.ix[i,'conc (ug/m^3)']
    elif C <1.5:
        histlist[3] += newarray.ix[i,'conc (ug/m^3)']
    elif C <2.5:
        histlist[4] += newarray.ix[i,'conc (ug/m^3)']
    elif C <3.5:
        histlist[5] += newarray.ix[i,'conc (ug/m^3)']
    elif C <4.5:
        histlist[6] += newarray.ix[i,'conc (ug/m^3)']
    elif C <5.5:
        histlist[7] += newarray.ix[i,'conc (ug/m^3)']
    elif C <20:
        histlist[8] += newarray.ix[i,'conc (ug/m^3)']      
        
with open(keyword+'histdata.txt','w') as a:
    for i in histlist:
        a.write(str(i)+'\n')
        
         
if os.path.isfile(root+'_aer_histdata.txt') and os.path.isfile(root+'_histdata.txt'):
    g=[]
    ae=[]   
    h1 = open(root+'_histdata.txt','r')
    for line in h1:
        #print line
        g.append(float(line.strip()))
    h1.close()
    h2 = open(root+'_aer_histdata.txt','r')
    for line in h2:
        ae.append(float(line.strip()))  
    h2.close()

    hist=pd.DataFrame(index=np.arange(-3,7,1),columns=['aer','gas'])
    hist['aer'] = ae
    hist['gas'] = g
        
    plot = hist.plot(title=root,legend=True,kind='bar',stacked=True,
                     color='gw',log=True)
    plt.legend(loc=2)
    plt.xlabel('log C^0')
    plt.xticks(rotation=0)
    plt.ylabel('ug/m^3')
    plt.ylim([0, 150])
    
    plt.show()

    plt.savefig(root+'satconcdist.pdf')


