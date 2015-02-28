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
parser.add_argument('mode',type=str,default='total') # {total|gas}

class InitModify:

    def __init__(self,root):
        self.root = root

        self.declaration = '''
    character(len=100) ::arg_string ! FB: Command-line arguments as string
    integer            :: narg      ! FB: # of args
    integer            :: keyboard_input
    integer            :: filestat
    logical            :: existp
    REAL(kind=dp)      :: DURATION  ! FB
    REAL(kind=dp)      :: CFACTOR_NEW=1.D0
    REAL(kind=dp)      :: CFACTOR_RATIO=1.D0
'''

        self.declaration_aer = '''
    USE apinene_GlobalAER, only: integratorcheck, partition_substeps, &
         cAER0_total, runpath ! FB
'''        
        
        self.define_inputs = '''
    !~~~> Parse command-line-arguments         ! FB
    !Check if any arguments are found          ! FB
    narg=command_argument_count()              ! FB
    if(narg>0)then                             ! FB
       call get_command_argument(1,arg_string) ! FB: Read in first argument (! Zeroth argument corresponds to program name) ! FB
       read(arg_string,"(I10)") keyboard_input ! FB
    else                                       ! FB
       keyboard_input = -1                     ! FB input_file number will be asked within the program (apinene_Initialize.f90) ! FB
    end if                                     ! FB ! FB
    !
    if (keyboard_input == -1) then
       ! FB decide whether to ovwerwrite initialization with data from external file
       print*,"Do you want to overwrite simulation parameters with data from an external file?"
       print*,"  - enter 999 to keep simulation parameters from ROOT.def and run_def/input.txt"
       print*,"  - or enter number of runpath (max 3 digits, run_DDD.txt) to overwrite definition from ROOT.def"
       read (*,*) keyboard_input
    end if

! define input.txt formats
110 format (f9.0)
111 format (i9)
112 format (a8)
113 format (E10.0)

    select case (keyboard_input)
    case(999) ! default
       runpath = "run_def"
    case DEFAULT ! FB
       write(runpath,"(A4,I0.3)") "run_",keyboard_input
    end select
'''

        self.read_aer = '''        
    ! partitioning (required)
    write(*,*) "using parameters from ", runpath//"/input_partitioning.txt"
    open (unit=15, file=runpath//"/input_partition.txt", status='old',    &
         access="sequential", form="formatted", action="read" )

    ! read in defaults defined in input.txt
    read (15, 113)  cAER0_total
    read (15, 111)  integratorcheck
    read (15, 111)  partition_substeps !0=no partitioning
    
    close(15)
    ! echo
    write(*,*) "read from file: ", runpath//"/input_partitioning.txt"
    write(*,*) "read in values: ", &
         "total initial CAER [g/m3]: ", cAER0_total, &
         "partition_substeps: ", partition_substeps, &
         "integratorcheck: ", integratorcheck 
'''
        self.read_optional = '''
    ! overwrite time (optional) (tested)
    inquire(file=runpath//"/input_time.txt", exist=existp)
    if (existp) then
       open (unit=15, file=runpath//"/input_time.txt", status='old',    &
             access="sequential", form="formatted", action="read" )

       read (15, 110)  TSTART
       read (15, 110)  DURATION
       read (15, 110)  DT

       close(15)

       TEND = TSTART + DURATION

       write(*,*) "read from file: ", runpath//"/input_time.txt"
       write(*,*) "read in: ", &
            "TSTART: ", TSTART, &
            "DURATION: ", DURATION, &
            "DT: ", DT 
    endif
        
    ! overwrite initial concentrations (optional) (untested)
    inquire(file=runpath//"/cgas_init.txt", exist=existp)
    if (existp) then
       open (unit=15, file=runpath//"/cgas_init.txt", status='old',    &
             access="sequential", form="formatted", action="read" )
       x = (1e-5)*CFACTOR
       do
          read (15, *, iostat=filestat)  idx, conc ! CHECK PRECISION
          if (filestat /= 0) exit
          VAR(idx) = conc*x
       end do
       close(15)

       TEND = TSTART + DURATION

       write(*,*) "read from file: ", runpath//"/input_time.txt"
       write(*,*) "read in: ", &
            "TSTART: ", TSTART, &
            "DURATION: ", DURATION, &
            "DT: ", DT 
    endif

    ! overwrite temperature (optional) (untested)
    inquire(file=runpath//"/input_temp.txt", exist=existp)
    if (existp) then
       open (unit=15, file=runpath//"/input_temp.txt", status='old',    &
             access="sequential", form="formatted", action="read" )

       read (15, 110)  TEMP
       read (15, 110)  CFACTOR_NEW ! make sure double precision?

       close(15)

       CFACTOR_RATIO = CFACTOR_NEW/CFACTOR
       CFACTOR = CFACTOR_NEW

       DO i = 1, NVAR
          VAR(i) = VAR(i)*CFACTOR_RATIO
       END DO

       DO i = 1, NFIX
          FIX(i) = FIX(i)*CFACTOR_RATIO
       END DO

       write(*,*) "read from file: ", runpath//"/input_temp.txt"
       write(*,*) "read in: ", &
            "CFACTOR: ", TSTART, &
            "DURATION: ", DURATION, &
            "DT: ", DT 
    endif
    ! End overwriting initialization FB/ST
'''.format(ROOT=self.root)

    def modify(self),mode:
        mode = True if mode=='total' else False
        filename = self.root+'_Initialize.f90'
        newfile = filename + '~'
        with open(newfile,'w') as fout:
            with open(filename) as finp:
                ## import modules and initialize
                for line in finp:
                    fout.write(line)
                    if 'USE' in line and '_Util' in line:
                        if mode:
                            fout.write(self.declaration_aer)
                        fout.write(self.declaration)
                        break
                ## overwrite INLINED initializations (i.e. leave them in the code but overwrite their result just after)
                for line in finp:
                    fout.write(line)
                    if 'End INLINED initializations' in line:
                        fout.write(self.define_inputs)
                        if mode:
                            fout.write(self.read_aer)
                        fout.write(self.read_optional)
                        fout.write('\n')
        os.rename(newfile,filename)
        print filename+' modified'

if __name__ == '__main__':

    args = parser.parse_args()
    root, mode = args.root, args.mode
    ## root = 'octane_gen'
    mainfile = InitModify(root, mode)
    mainfile.modify()
