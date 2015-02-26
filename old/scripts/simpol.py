######################################################################
##* Script to calculate group abundances and vapor pressures
##* by SIMPOL. Group abundances enumerated by matching of
##* SMILES and SMARTS patterns using pybel.
##* ~simpol.py~
##* $Rev$
##* Sep. 2013
##* Satoshi Takahama (satoshi.takahama@epfl.ch)
######################################################################

import pandas as pd
import numpy as np
import cStringIO
import pybel
from collections import OrderedDict

def mergedict(self,*args):
    out = {}
    for d in args:
        out.update(d)
    return out

def desblet(args,expr):
    ## combines destructuring bind with let
    locals().update(args)
    return eval(expr)

class Simpolclass:
    def __init__(self):
        # from Pankow and Asher (2008)
        TABLE5 = """groups	k	coefficient	footnote comment	Bk,1	Bk,2	Bk,3	Bk,4
zeroeth group (constant term)	0	b0	a	-4.26938E+02	2.89223E-01	4.42057E-03	2.92846E-01
carbon number	1	b1	b	-4.11248E+02	8.96919E-01	-2.48607E-03	1.40312E-01
carbon number on the acid-side of an amide (asa)	2	b2	c	-1.46442E+02	1.54528E+00	1.71021E-03	-2.78291E-01
aromatic ring	3	b3	d	3.50262E+01	-9.20839E-01	2.24399E-03	-9.36300E-02
non-aromatic ring	4	b4	e	-8.72770E+01	1.78059E+00	-3.07187E-03	-1.04341E-01
C=C (non-aromatic)	5	b5	f	5.73335E+00	1.69764E-02	-6.28957E-04	7.55434E-03
C=C-C=O in non-aromatic ring	6	b6	g	-2.61268E+02	-7.63282E-01	-1.68213E-03	2.89038E-01
hydroxyl (alkyl)	7	b7	h	-7.25373E+02	8.26326E-01	2.50957E-03	-2.32304E-01
aldehyde	8	b8	i	-7.29501E+02	9.86017E-01	-2.92664E-03	1.78077E-01
ketone	9	b9	j	-1.37456E+01	5.23486E-01	5.50298E-04	-2.76950E-01
carboxylic acid	10	b10	k	-7.98796E+02	-1.09436E+00	5.24132E-03	-2.28040E-01
ester	11	b11	L	-3.93345E+02	-9.51778E-01	-2.19071E-03	3.05843E-01
ether	12	b12	m	-1.44334E+02	-1.85617E+00	-2.37491E-05	2.88290E-01
ether (alicyclic)	13	b13	m	4.05265E+01	-2.43780E+00	3.60133E-03	9.86422E-02
ether, aromatic	14	b14	m	-7.07406E+01	-1.06674E+00	3.73104E-03	-1.44003E-01
nitrate	15	b15	n	-7.83648E+02	-1.03439E+00	-1.07148E-03	3.15535E-01
nitro	16	b16	o	-5.63872E+02	-7.18416E-01	2.63016E-03	-4.99470E-02
aromatic hydroxyl (e.g., phenol)	17	b17	p	-4.53961E+02	-3.26105E-01	-1.39780E-04	-3.93916E-02
amine, primary	18	b18	q	3.71375E+01	-2.66753E+00	1.01483E-03	2.14233E-01
amine, secondary	19	b19	q	-5.03710E+02	1.04092E+00	-4.12746E-03	1.82790E-01
amine, tertiary	20	b20	q	-3.59763E+01	-4.08458E-01	1.67264E-03	-9.98919E-02
amine, aromatic	21	b21	q	-6.09432E+02	1.50436E+00	-9.09024E-04	-1.35495E-01
amide, primary	22	b22	c	-1.02367E+02	-7.16253E-01	-2.90670E-04	-5.88556E-01
amide, secondary	23	b23	c	-1.93802E+03	6.48262E-01	1.73245E-03	3.47940E-02
amide, tertiary	24	b24	c	-5.26919E+00	3.06435E-01	3.25397E-03	-6.81506E-01
carbonylperoxynitrate	25	b25	r	-2.84042E+02	-6.25424E-01	-8.22474E-04	-8.80549E-02
peroxide	26	b26	r	1.50093E+02	2.39875E-02	-3.37969E-03	1.52789E-02
hydroperoxide	27	b27	r	-2.03387E+01	-5.48718E+00	8.39075E-03	1.07884E-01
carbonylperoxyacid	28	b28	r	-8.38064E+02	-1.09600E+00	-4.24385E-04	2.81812E-01
nitrophenol	29	b29	p	-5.27934E+01	-4.63689E-01	-5.11647E-03	3.84965E-01
nitroester	30	b30	L	-1.61520E+03	9.01669E-01	1.44536E-03	2.66889E-01
"""
        self.table = pd.read_table(cStringIO.StringIO(TABLE5))
        ## SMARTS patterns required by SIMPOL
        self.smartspatt = OrderedDict([
            ('_zero','_zero'),
            ('_carbon number','[#6]'),
            ('_carbon number on the acid-side of an amide (asa)','_nomatch'),
            ('_aromatic ring','_nomatch'),            
            ('_non-aromatic ring','_nomatch'),
            ('_C=C (non-aromatic)','C=C'),
            ('_C=C-C=O in non-aromatic ring','[$(C=CC=O);A;R]'),
            ('_hydroxyl (alkyl)','[C;!$(C=O)][OX2H,OD1]'),#includes radicals
            ('_aldehyde','[CX3H1](=O)[#6]'),
            ('_ketone','[#6][CX3](=O)[#6]'),
            ('_carboxylic acid','[CX3](=O)[OX2H,OD1]'),#includes radicals
            ('_ester','[CX3,CX3H1](=O)[OX2H0][#6]'),
            ('_ether','_nomatch'),
            ('_ether (alicyclic)','[OD2;R]([C;!$(C=O)])[C;!$(C=O)]'),
            ('_ether, aromatic','c[OD2][c,C&!$(C=O)]'),
            ('_nitrate','[$([NX3](=[OX1])(=[OX1])[O;$([X2]),$([X1-])]),$([NX3+]([OX1-])(=[OX1])[O;$([X2]),$([X1-])])]'),
            ('_nitro','[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'),
            ('_aromatic hydroxyl','c[OX2H,OD1]'),#includes radicals
            ('_amine, primary','[NX3;H2;!$(NC=O);!$(Na)]'),
            ('_amine, secondary','[NX3;H1;!$(NC=O);!$(Na)]'),
            ('_amine, tertiary','[NX3;H0;!$(NC=O);!$(N=O);!$(Na)]'),
            ('_amine, aromatic','[N;!$(NC=O);!$(N=O);$(Na)]'),
            ('_amide, primary','[CX3](=O)[NX3H2]'),
            ('_amide, secondary','[CX3](=O)[NX3H1]'),
            ('_amide, tertiary','[CX3](=O)[NX3H0;!$(N=O)]'),
            ('_carbonylperoxynitrate','C(=O)OO[$([NX3](=[OX1])(=[OX1])[O;$([X2]),$([X1-])]),$([NX3+]([OX1-])(=[OX1])[O;$([X2]),$([X1-])])]'),
            ('_peroxide','[OD2][OD2,OD1]'),
            ('_hydroperoxide','[OX2H,OD1][OD2][#6;!$(C=O)]'),#includes radicals
            ('_carbonylperoxyacid','C(=O)O[OX2H,OD1]'),#includes radicals
            ('_nitrophenol','[O-][N+](=O)[$(cc(O)cccc),$(ccc(O)ccc),$(cccc(O)cc)]'),
            ('_nitroester','[#6][OX2H0][CX3,CX3H1](=O)[C;$(C[N+](=O)[O-]),$(CC[N+](=O)[O-]),$(CCC[N+](=O)[O-]),$(CCCC[N+](=O)[O-]),$(CCCCC[N+](=O)[O-])]')
           # ('_alkane CH', '_nomatch'),            
           # ('_tertiary C','[CX4H1]'),
            ])
        self.smartsaux = OrderedDict([
            ('_ac5','[a;r5]'),
            ('_ac6','[a;r6]'),
            ('_AC3','[A;r3]'),
            ('_AC4','[A;r4]'),
            ('_AC5','[A;r5]'),
            ('_AC6','[A;r6]'),
            ('_AC7','[A;r7]'),
            ('_AC8','[A;r8]'),
            ('_AC9','[A;r9]'),
            ('_AC10','[A;r10]'),
            ('_ether_all','[OD2]([#6])[#6]'),
            ('_asa aux','[C;$(C[NX3][CH,CC](=O)),$(CC[NX3][CH,CC](=O)),$(CCC[NX3][CH,CC](=O)),$(CCCC[NX3][CH,CC](=O)),$(CCCCC[NX3][CH,CC](=O))]'),
            ('_CH2','[CX4H2]'),
            ('_CH3','[CX4H3]'),
            ('_edges','*~*'),
            ('_atoms','*')
            ])
        self.smartsmath = OrderedDict([
            ('','')
            ])


    def p0_atm(self,nuk,temp):
        ## calculates vapor pressure from vector of abundances (nuk) and
        ##   temperature (temp in K)
        bk = self.table['Bk,1']/temp + self.table['Bk,2']+\
             self.table['Bk,3']*temp + self.table['Bk,4']*np.log(temp)
        return 10**np.sum(nuk*bk)

    def deltaHvap_kJpermol(self,nuk,temp):
        ## calculates enthalpy of vaporization from vector of abundances (nuk) and
        ##   temperature (temp in K)
        ## R is the gas constant
        ## db is db(T)/d(1/T)
        R = 8.3144622e-3 ## kJ K^-1 mol^-1
        dbT = self.table['Bk,1']-self.table['Bk,3']*temp**2-self.table['Bk,4']*temp
        return -2.303 * R * np.sum(nuk*dbT)

    def press2conc(self,press,temp):
        ## converts pressure (atm) to concentration (molecules/cm^3)
        Av = 6.0221413e+23 # molecules/mole
        R = 82.05746       # cm^3 atm/ mole/ K
        conc = press/temp * Av/R
        return conc

    def get_groupnames(self):
	return self.table['groups'].tolist()

    def get_groups(self,smilesstr=None):
        ## takes a single SMILES string and returns an array (length=31) 
        ##   of group abundances for SIMPOL
        ##
        ## if smilesstr is none, return array of zeroes
        if not smilesstr:
            return np.array([(1 if '_zero'==p else 0) for p in self.smartspatt.keys()])
        ## main body
        mol = pybel.readstring('smi',smilesstr)
        ## store SIMPOL patterns in ordered dictionary
        abundances = OrderedDict()
        abundances['_zero'] = 1
        for key,patt in self.smartspatt.items()[1:]:
            abundances[key] = 0 if 'nomatch' in patt else \
                              len(pybel.Smarts(patt).findall(mol))
        ## find auxiliary patterns
        aux = OrderedDict()
        for key,patt in self.smartsaux.items():
            aux[key] = len(pybel.Smarts(patt).findall(mol))
        ## SMARTS math
        combineddict = mergedict(abundances,aux)
        for key,patt in self.smartsmath.items():
            if patt=='':
                continue
            abundances[key] = desblet(combineddict,patt)
        
	#nonaromatic rings (AC)
        Aringcount = 0
        for ringsize in range(3,11): #for ringsizes from 3 to 10
            key = '_AC'+str(ringsize)
            Aringcount += aux[key]/ringsize + (aux[key] % ringsize > 0)
        abundances['_non-aromatic ring'] = Aringcount


        abundances['_aromatic ring'] = aux['_edges'] - aux['_atoms'] + 1 - abundances['_non-aromatic ring']

        
        abundances['_ether'] = aux['_ether_all'] \
                             - abundances['_ester'] \
                             - abundances['_ether (alicyclic)'] \
                             - abundances['_ether, aromatic']
                             
        abundances['_nitrate'] = abundances['_nitrate'] \
                               - abundances['_carbonylperoxynitrate']
                               
        #assumes only 1 nitrophenol ring per molecule and no nitrodiphenols
        if abundances['_nitrophenol'] > 0:
            abundances['_nitrophenol'] = 1                               
                               
        abundances['_aromatic hydroxyl'] = abundances['_aromatic hydroxyl'] \
                                         - abundances['_nitrophenol']
                                         
        abundances['_ester'] -= abundances['_nitroester']                                         

        abundances['_peroxide'] = abundances['_peroxide'] \
                                - abundances['_carbonylperoxynitrate'] \
                                - abundances['_hydroperoxide'] \
                                - abundances['_carbonylperoxyacid']
       
        if abundances['_amide, primary'] > 0 or \
           abundances['_amide, secondary'] > 0 or \
           abundances['_amide, tertiary'] > 0:
            abundances['_carbon number on the acid-side of an amide (asa)'] = \
            abundances['_carbon number'] - aux['_asa aux'] - 1
        
       # abundances['_alkane CH'] = abundances['_tertiary C'] + aux['_CH2']*2 + aux['_CH3']*3
        
	return np.array(abundances.values())
    
    def calc_properties(self,smiles,temp):
        ## calculate vapor pressures and Delta H for list of smiles strings
        ## return value is a Series obj
        if type(smiles) is str: # if single string (not list)
            smiles = [smiles]
        props = pd.DataFrame(np.tile(np.nan,(len(smiles),2)),
                             index=smiles.index,
                             columns=['p0','DeltaH'])
        for name in smiles.index:
            abundances = self.get_groups(smiles[name])
            if any(abundances[1:] > 0):
                props.ix[name,'p0'] = self.p0_atm(abundances,temp)
                props.ix[name,'DeltaH'] = self.deltaHvap_kJpermol(abundances,temp)
        return props

if __name__=='__main__':

    ## run example
    test = pd.Series(OrderedDict([
        ('NC8H18', 'CCCCCCCC'),
        ('CH3CHO', 'CC=O'),
        ('APINENE', 'CC1=CCC2CC1C2(C)C')
        ]))
    simp = Simpolclass()
    print simp.get_groups(test[0]) ## abundances
    print simp.get_groups(test[2]) ## abundances    
    print simp.calc_properties(test,298.15) ## vapor pressures (atm)
