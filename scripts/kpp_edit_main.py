#!/usr/bin/env python

################################################################################
##
## kpp_edit_main.py
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
################################################################################

import os
import argparse

parser = argparse.ArgumentParser(description='edit {root}_Main.f90 file')
parser.add_argument('root',type=str)

class MainDriver:

    def __init__(self,root):
        self.root = root
        self.importmod = '''  USE {ROOT}_Partition, ONLY: PARTITION !ST
  USE {ROOT}_InitializeAER, ONLY: InitializeAER                     !ST
  USE {ROOT}_UtilAER                                                ! ST
  use {ROOT}_GlobalAER, only: partition_substeps                    ! FB  
'''.format(ROOT=self.root)

    def modify(self):
        filename = self.root+'_Main.f90'
        newfile = filename + '~'
        addstatement = lambda x,n: ' '*n+x+' ! ST\n'
        with open(newfile,'w') as fout:
            with open(filename) as finp:
                ## import modules
                for line in finp:
                    fout.write(line)
                    if 'USE' in line and '_Initialize' in line:
                        fout.write(self.importmod)
                        break
                ## initialize
                for line in finp:
                    fout.write(line)
                    if 'CALL InitSaveData()' in line:
                        fout.write(addstatement('CALL InitializeAER()',6))
                        fout.write(addstatement('CALL InitSaveDataAER()',6))
                        break
                ## save data, call partition
                accum = ''
                for line in finp:
                    if 'CALL SaveData()' in line:
                        accum += line
                        accum += addstatement('CALL SaveDataAER()',8)
                        next
                    elif 'END DO kron' in line:
                        accum += addstatement('IF (PARTITION_SUBSTEPS > 0) THEN',8)
                        accum += addstatement('TNEXT=T+DT',10)
                        accum += addstatement('CALL PARTITION(T,TNEXT)',10)
                        accum += addstatement('ENDIF',8)
                        accum += line
                        break
                    accum += line
                fout.write(accum)
                ## save data
                for line in finp:
                    fout.write(line)
                    if 'CALL CloseSaveData()' in line:
                        fout.write(addstatement('CALL SaveDataAER()',6))
                        fout.write(addstatement('CALL CloseSaveDataAER()',6))
        os.rename(newfile,filename)
        print filename+' modified'

if __name__ == '__main__':

    args = parser.parse_args()
    root = args.root
    ## root = 'octane_gen'
    mainfile = MainDriver(root)
    mainfile.modify()
