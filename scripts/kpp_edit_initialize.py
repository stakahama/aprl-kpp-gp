#!/usr/bin/env python

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

        self.module_aer = '''
    USE apinene_GlobalAER, only: cAER0_total, &
         integratorcheck, partition_on, minconc, mf
'''

        self.declaration_gas = '''
    REAL(kind=dp)      :: minconc             ! FB
'''

        
        self.declaration = '''
    integer            :: filestat            ! ST
    logical            :: existp              ! ST
    REAL(kind=dp)      :: DURATION            ! FB
    REAL(kind=dp)      :: CFACTOR_NEW=1.D0    ! ST
    REAL(kind=dp)      :: CFACTOR_RATIO=1.D0  ! ST
    integer            :: ix                  ! ST
    REAL(kind=dp)      :: conc                ! ST    
'''

        self.define_inputs = '''
        
! define input.txt formats for inputs
110 format (f9.0)
111 format (i9)
112 format (a8)
113 format (E10.0)

'''

        self.minconc = '''
    minconc = 1.0D-5*CFACTOR
'''

        self.read_aer = '''        
    ! partitioning (required)
    write(*,*) "using parameters from ", "input_partitioning.txt"
    open (unit=15, file="input_partitioning.txt", status="old",    &
         access="sequential", form="formatted", action="read")

    ! read in defaults defined in input.txt
    read (15, 113)  cAER0_total        !initial aerosol concentration
    read (15, 111)  partition_on       !0=off; 1=on
    read (15, 111)  integratorcheck    !0=off; 1=on
    read (15, 113)  minconc            !minimum concentration (zero value) [ppb]
    read (15, 111)  mf                 !LSODE Jacobian option (10, 21, 22)
    
    close(15)
    ! echo
    write(*,*) "read from file: ", "input_partitioning.txt"
    write(*,*) "read in values: ", &
         "total initial CAER [g/m3]: ", cAER0_total, &
         "partition_on: ", partition_on, &
         "integratorcheck: ", integratorcheck, &
         "minconc: ", minconc, &          
         "mf: ", mf         

    minconc = minconc*CFACTOR         
'''
        self.read_optional = '''
    ! overwrite time (optional)
    inquire(file="input_time.txt", exist=existp)
    if (existp) then
       open (unit=15, file="input_time.txt", status='old',    &
             access="sequential", form="formatted", action="read")

       read (15, 110)  TSTART
       read (15, 110)  DURATION
       read (15, 110)  DT

       close(15)

       TEND = TSTART + DURATION

       write(*,*) "read from file: ", "input_time.txt"
       write(*,*) "read in: ", &
            "TSTART: ", TSTART, &
            "DURATION: ", DURATION, &
            "DT: ", DT 
    endif
        
    ! overwrite initial concentrations
    inquire(file="cgas_init.txt", exist=existp)
    if (existp) then
       open (unit=15, file="cgas_init.txt", status="old",    &
             access="sequential", form="formatted", action="read")
       do
          read (15, *, iostat=filestat)  ix, conc ! CHECK PRECISION
          if (filestat /= 0) exit
          if (ix .eq. 0) then
             H2O = conc
          else       
             VAR(ix) = conc*CFACTOR
          endif
       end do
       close(15)

       write(*,*) "read from file: ", "cgas_init.txt"
    endif

    ! overwrite temperature (optional) (untested)
    inquire(file="input_temp.txt", exist=existp)
    if (existp) then
       open (unit=15, file="input_temp.txt", status="old",    &
             access="sequential", form="formatted", action="read")

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

       minconc = minconc*CFACTOR_RATIO

       write(*,*) "read from file: ", "input_temp.txt"
       write(*,*) "read in: ", &
            "TEMP (K): ", TEMP, &
            "CFACTOR (molec/cm^3/ppb): ", CFACTOR_NEW
    endif
    ! End overwriting initialization FB/ST
'''.format(ROOT=self.root)

    def modify(self,mode):
        modeTF = {'gas':False,'total':True}
        mode = modeTF[mode]
        filename = self.root+'_Initialize.f90'
        newfile = filename + '~'
        with open(newfile,'w') as fout:
            with open(filename) as finp:
                ## import modules and initialize
                for line in finp:
                    fout.write(line)
                    if 'USE' in line and '_Util' in line:
                        if mode:
                            fout.write(self.module_aer)
                        else:
                            fout.write(self.declaration_gas)
                        fout.write(self.declaration)
                        break
                ## overwrite INLINED initializations (i.e. leave them in the code but overwrite their result just after)
                for line in finp:
                    fout.write(line)
                    if 'End INLINED initializations' in line:
                        fout.write(self.define_inputs)
                        fout.write(self.minconc)                        
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
    mainfile = InitModify(root)
    mainfile.modify(mode)
