#!/usr/bin/env python

################################################################################
##
## kpp_edit_makefile.py
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
## license: GNU Public License v3.0 (LICENSE_GPLv3.txt)
##
################################################################################

import os
import argparse

parser = argparse.ArgumentParser(description='edit Makefile_{root}.f90')
parser.add_argument('root',type=str)

class NewMakeDef:

    def __init__(self, root):

        self.root = root

        self.aerobj = r'''## ST
AEROBJ = simpol_module.o \
         {ROOT}_GlobalAER.o \
	 {ROOT}_InitializeAER.o \
	 {ROOT}_UtilAER.o \
	 {ROOT}_PartitionAER.o \
	 {ROOT}_SIMPOLGroups.o \
         opkdmain.o \
         opkda1.o \
         opkda2.o
## end ST
'''.format(ROOT=self.root)

        self.dotofiles = '''
## ST
{ROOT}_GlobalAER.o: {ROOT}_GlobalAER.f90 $(GENOBJ)
	$(FC) $(FOPT) -x f95-cpp-input -c $<

{ROOT}_InitializeAER.o: {ROOT}_InitializeAER.f90 $(GENOBJ) {ROOT}_GlobalAER.o {ROOT}_SIMPOLGroups.o
	$(FC) $(FOPT) -x f95-cpp-input -c $<

{ROOT}_UtilAER.o: {ROOT}_UtilAER.f90 $(GENOBJ) $(UTLOBJ) {ROOT}_GlobalAER.o
	$(FC) $(FOPT) -x f95-cpp-input -c $<

{ROOT}_PartitionAER.o: {ROOT}_PartitionAER.f90 $(GENOBJ) {ROOT}_GlobalAER.o
	$(FC) $(FOPT) -x f95-cpp-input -c $<

{ROOT}_SIMPOLGroups.o: {ROOT}_SIMPOLGroups.f90
	$(FC) $(FOPT) -x f95-cpp-input -c $<

constants.o: kpp_constants.f90
	$(FC) $(FOPT) -x f95-cpp-input -c $< -o $@

simpol_module.o: simpol_module.f90
	$(FC) $(FOPT) -x f95-cpp-input -c $<

opkdmain.o: opkdmain.f
	$(FC) $(FOPT) -x f77-cpp-input -c $<

opkda1.o: opkda1.f
	$(FC) -x f77-cpp-input -c $<

opkda2.o: opkda2.f
	$(FC) -x f77-cpp-input -c $<
## end ST
'''.format(ROOT=self.root)

    def modify(self):
        ##
        filename = filename = 'Makefile_' + self.root
        newfile = filename+'~'
	eol = '\n'
        ##
        append2ALLOBJ = lambda x: x.replace(eol,' $(AEROBJ)'+eol)
        append2exe = lambda x: x.replace(' -o',' constants.o -o')
	append2clean = lambda x: x.replace('\\','constants.o constants.mod simpol_module.o simpol_module.mod *histdata.txt \\')
        replrates = lambda x: x.replace(eol,' constants.o'+eol)
        ##
        with open(newfile,'w') as fout:
            with open(filename) as finp:
            ## define aerobj
                for line in finp:
                    fout.write(line)
                    if 'GENOBJ =' == line[:8]:
                        break
                for line in finp:
                    fout.write(line)
                    if line.strip()=='':
                        break
                fout.write(self.aerobj+'\n')
            ## add aerobj
                for line in finp:
                    fout.write(line)
                    if 'ALLOBJ =' == line[:8]:
                        break
                fout.write(append2ALLOBJ(next(finp)))
            ## add kpp_constants.o
                for line in finp:
                    fout.write(line)
                    if 'exe:' == line[:4]:
                        break
                fout.write(append2exe(next(finp)))
	    ## add to clean and distclean: constants.mod, simpol_module.mod
                for line in finp:
                    fout.write(line)
                    if 'clean:' == line[:6]:
                        break
                fout.write(append2clean(next(finp)))
                for line in finp:
                    fout.write(line)
                    if 'distclean:' == line[:10]:
                        break
                fout.write(append2clean(next(finp)))
            ## add constants.o to rates
                for line in finp:
                    if '_Rates.o:' not in line:
                        fout.write(line)
                    else:
                        fout.write(replrates(line))
                        break
            ## write out rest
                for line in finp:
                    fout.write(line)
            ## define targets
            fout.write(self.dotofiles+'\n')
        ## overwrite
        os.rename(newfile,filename)
        print filename+' modified'

if __name__=='__main__':

    args = parser.parse_args()
    root = args.root
    ## root = 'octane_gen'
    makefile = NewMakeDef(root)
    makefile.modify()
