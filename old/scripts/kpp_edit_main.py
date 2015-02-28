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
  USE {ROOT}_InitializeAER, ONLY: InitializeAER !ST
  USE {ROOT}_UtilAER ! ST
  use {ROOT}_Global, only: TEMP, NSPEC, C              ! FB
  use {ROOT}_GlobalAER, only: integrator, partition_substeps, keyboard_input   ! FB
  character(len=100)::arg_string ! FB: Command-line arguments as string
  integer::narg ! FB: # of args
'''.format(ROOT=self.root)

    def modify(self):
        filename = self.root+'_Main.f90'
        newfile = filename + '~'
        addstatement = lambda x,n: ' '*n+x+' ! ST\n'
        addstatementFB = lambda x,n: ' '*n+x+' ! FB\n'
        with open(newfile,'w') as fout:
            with open(filename) as finp:
                ## import modules
                for line in finp:
                    fout.write(line)
                    if 'USE' in line and '_Initialize' in line:
                        fout.write(self.importmod)
                        break
                ## parse command line arguments
                for line in finp:
                    fout.write(line)
                    if 'INTEGER :: i' in line:
                        fout.write(addstatementFB('!~~~> Parse command-line-arguments          ',0))
                        fout.write(addstatementFB('  !Check if any arguments are found         ',0))
                        fout.write(addstatementFB('   narg=command_argument_count()            ',0))
                        fout.write(addstatementFB('   if(narg>0)then                           ',0))
                        fout.write(addstatementFB('    call get_command_argument(1,arg_string) ! FB: Read in first argument (! Zeroth argument corresponds to program name)',0))
                        fout.write(addstatementFB('    read(arg_string,\'(I10)\') keyboard_input ',0))                            
                        fout.write(addstatementFB('  else                                      ',0))                            
                        fout.write(addstatementFB('    keyboard_input = -1  ! FB input_file number will be asked within the program ({ROOT}_Initialize.f90)'.format(ROOT=self.root),0))                            
                        fout.write(addstatementFB('  end if                                    ! FB',0))                            
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
                        continue
                    elif 'END DO kron' in line:
                        #accum += addstatement('CALL PARTITION()',8)+'\n'
                        accum += addstatementFB('select case (integrator)',8)
                        accum += addstatementFB('  case ("backward","forward","dlsode")',8) # TODO (FB, 16.12.14): check whether these cases still work, I replaced the ' with " when writing it into the python script
                        accum += addstatementFB('    CALL PARTITION(DT, partition_substeps) ! ST & FB 1 = number of intermediate (sub-)timesteps to interpolate',8)
                        accum += addstatementFB('  case DEFAULT',8)
                        accum += addstatementFB('    write(*,*) "no partitioning to aerosol selected"',8)
                        accum += addstatementFB('    CALL PARTITION(DT, 0) ! ST & FB 1 = number of intermediate (sub-)timesteps to interpolate, NOTE: I still call the function in order to have the output',8)
                        accum += addstatementFB('end select',8)
                        accum += '\n'
                        accum += line
                        break
                    accum += line
                fout.write(accum)
                ## save data
                for line in finp:
                    fout.write(line)
                    if 'CALL CloseSaveData()' in line:
                        # TODO (FB, 16.12.14): remove, just for printing to inputXXoutput_SCREEN.txt
                        accum += addstatementFB('WRITE(88,991) (T-TSTART)/(TEND-TSTART)*100, T,     &',8)
                        accum += addstatementFB('         ( TRIM(SPC_NAMES(MONITOR(i))),           &',8)
                        accum += addstatementFB('           C(MONITOR(i))/CFACTOR, i=1,NMONITOR ), &',8)
                        accum += addstatementFB('         ( TRIM(SMASS(i)), DVAL(i)/CFACTOR, i=1,NMASS )',8)
                        # END TODO
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
