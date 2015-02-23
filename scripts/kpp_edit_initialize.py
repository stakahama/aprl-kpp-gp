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
        self.init_modification = '''if (keyboard_input == -1) then
! FB decide whether to ovwerwrite initialization with data from external file
print*,"Do you want to overwrite simulation parameters with data from an external file?"
print*,"  - enter 999 to keep simulation parameters from ROOT.def and input_template_.txt"
print*,"  - or enter number of input file (max 2 digits, inputXX.txt) to overwrite definition from ROOT.def"
read (*,*) keyboard_input
end if

! define input_XX.txt formats
110 format (f9.0)
111 format (i9)
112 format (a8)
113 format (E10.0)

select case (keyboard_input)
  case(999)
    write(*,*) "using parameters from ROOT.def and input_template_.txt"
    open (unit=15, file="in_outputs/input_template_.txt", status='old',    &
                 access='sequential', form='formatted', action='read' )
      ! read in defaults defined in input_template_.txt
        read (15, *)            ! skip first four lines
        read (15, *)  
        read (15, *)  
        read (15, *)  
        read (15, 111)  partition_substeps
        read (15, 112)  integrator
        read (15, 113)  cAER0_total

    input_file = "std-def"
    write(*,*) "read from file: ","input_template.txt", ", write output to: ", input_file
    write(*,*) "using: TSTART: ", TSTART, "TEND: ", TEND, "DT: ", DT, &
     "TEMP: ", TEMP, "read in: partition_substeps: ", partition_substeps, "integrator: ", &
     integrator, "total initial CAER [g/m3]", cAER0_total

  case DEFAULT ! FB
    write(input_file,"(A5,I0.2)") "input",keyboard_input
    open (unit=15, file="in_outputs/"//input_file//".txt", status='old',    &
                 access='sequential', form='formatted', action='read' )
        
        read (15, 110)  TSTART
        read (15, 110)  DURATION
        read (15, 110)  DT
        read (15, 110)  TEMP
        read (15, 111)  partition_substeps
        read (15, 112)  integrator
        read (15, 113)  cAER0_total

    TEND = TSTART + DURATION

    write(*,*) "read from file: ", input_file
    write(*,*) "read in: TSTART: ", TSTART, "DURATION: ", DURATION, "DT: ", DT, &
     "TEMP: ", TEMP, "partition_substeps: ", partition_substeps, "integrator: ", &
     integrator, "total initial CAER [g/m3]", cAER0_total
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
                        fout.write(addstatementFB('USE {ROOT}_GlobalAER, only: integrator, partition_substeps, keyboard_input, input_file, cAER0_total'.format(ROOT=self.root),2))
                        fout.write(addstatementFB('REAL(kind=dp) :: DURATION',2))
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
