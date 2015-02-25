################################################################################
##
## kpp_edit_makefile.py
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
################################################################################

import os
import argparse

parser = argparse.ArgumentParser(description='edit Makefile_{root}.f90')
parser.add_argument('root',type=str)

class NewMakeDef:

    def __init__(self, root):

        self.root = root

        self.dotofiles = '''
## ST
constants.o: kpp_constants.f90
	$(FC) $(FOPT) -x f95-cpp-input -c $< -o $@
## end ST
'''.format(ROOT=self.root)

    def modify(self):
        ##
        filename = filename = 'Makefile_' + self.root
        newfile = filename+'~'
	eol = '\n'
        ##
        append2exe = lambda x: x.replace(' -o',' constants.o -o')
	append2clean = lambda x: x.replace('\\','constants.o constants.mod \\')
        replrates = lambda x: x.replace(eol,' constants.o'+eol)
        ##
        with open(newfile,'w') as fout:
            with open(filename) as finp:
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
