################################################################################
##
## kpp_edit_initialize
## F.Bernhard (fabian.bernhard@epfl.ch)
## December 2014
##
## adapted from: "kpp_edit_main.py"
## S. Takahama (satoshi.takahama@epfl.ch)
## June 2014
##
################################################################################

import os
import argparse

parser = argparse.ArgumentParser(description='edit {root}_Initialize.f90 file')
parser.add_argument('root',type=str)

class MainDriver:

    def __init__(self,root):
        self.root = root

        self.declarations = '''
  REAL(kind=dp)      :: DURATION             !
  character(len=100) :: arg_string           ! Command-line arguments as string  
  integer            :: keyboard_input, narg ! run number
  character(len=7)   :: runpath              ! run path
'''
        
        self.init_modification = '''
!~~~> Parse command-line-arguments           ! FB
!Check if any arguments are found            ! FB
   narg=command_argument_count()             ! FB
   if(narg>0)then                            ! FB
    call get_command_argument(1,arg_string) ! FB: Read in first argument (! Zeroth argument corresponds to program name) ! FB
    read(arg_string,'(I10)') keyboard_input  ! FB
  else                                       ! FB
    keyboard_input = -1  ! FB input_file number will be asked within the program (apinene_Initialize.f90) ! FB
  end if                                    ! FB ! FB
  
  if (keyboard_input == -1) then
  ! FB decide whether to ovwerwrite initialization with data from external file
  print*,"Do you want to overwrite simulation parameters with data from an external file?"
  print*,"  - enter 999 to keep simulation parameters from {ROOT}.def and run_def/input.txt"
  print*,"  - or enter number of runpath (max 2 digits, run_DD.txt) to overwrite definition from {ROOT}.def"
  read (*,*) keyboard_input
  end if

  ! define input.txt formats
110 format (f9.0)
111 format (i9)
112 format (a8)
113 format (E10.0)

  select case (keyboard_input)
    case(999)
      write(*,*) "using default parameters"

    case DEFAULT ! FB
      write(runpath,"(A4,I0.3)") "run_",keyboard_input
      open (unit=15, file=runpath//"/input.txt", status='old',    &
                   access='sequential', form='formatted', action='read' )
          
          read (15, 110)  TSTART
          read (15, 110)  DURATION
          read (15, 110)  DT
          read (15, 110)  TEMP

      TEND = TSTART + DURATION

      write(*,*) "read from file: ", runpath
      write(*,*) "read in: TSTART: ", TSTART, "DURATION: ", DURATION, "DT: ", DT, &
       "TEMP: ", TEMP
  end select
! End overwriting initialization FB
'''.format(ROOT=self.root)

    def modify(self):
        filename = self.root+'_Initialize.f90'
        newfile = filename + '~'
        addstatement = lambda x,n: ' '*n+x+' ! ST\n'
        addstatementFB = lambda x,n: ' '*n+x+' ! FB\n'
        with open(newfile,'w') as fout:
            with open(filename) as finp:
                ## import modules and initialize
                for line in finp:
                    fout.write(line)
                    if 'USE' in line and '_Util' in line:
                        fout.write(self.declarations)                        
                        fout.write('\n')
                        break
                ## overwrite INLINED initializations (i.e. leave them in the code but overwrite their result just after)
                for line in finp:
                    fout.write(line)
                    if 'End INLINED initializations' in line:
                        fout.write(self.init_modification)                        
                        fout.write('\n')
        os.rename(newfile,filename)
        print filename+' modified'

if __name__ == '__main__':

    args = parser.parse_args()
    root = args.root
    ## root = 'octane_gen'
    mainfile = MainDriver(root)
    mainfile.modify()
