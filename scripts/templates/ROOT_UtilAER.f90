module {ROOT}_UtilAER

  IMPLICIT NONE
  
  integer, parameter :: fileidAER = 34

contains

! ----- copied from {ROOT}_Utils -----

! ****************************************************************
!                            
! InitSaveDataAER - Open the data file and write header
!   Parameters :                                                  
!
! ****************************************************************

  SUBROUTINE InitSaveDataAER ()
! modules
    USE {ROOT}_Parameters
    USE {ROOT}_Monitor

! local variables
    INTEGER i

! body
    open(fileidAER, file='{ROOT}_aer.dat')

    WRITE(fileidAER,998) 'TIME',  &
         (TRIM(SPC_NAMES(LOOKAT(i))), i=1,NLOOKAT)
998 FORMAT(A25,100(1X,A25))

  END SUBROUTINE InitSaveDataAER

! ****************************************************************
!                            
! SaveDataAER - Open the data file and write header
!   Parameters :                                                  
!
! ****************************************************************

  SUBROUTINE SaveDataAER ()
! modules
    USE {ROOT}_Global
    USE {ROOT}_Monitor
    USE {ROOT}_Parameters
    USE {ROOT}_GlobalAER, only: CAER
! local variables
    INTEGER i
! body
    WRITE(fileidAER,999) (TIME-TSTART)/3600.D0,  &
         (CAER(LOOKAT(i))/CFACTOR, i=1,NLOOKAT)
999 FORMAT(ES25.16E3,100(1X,ES25.16E3))

  END SUBROUTINE SaveDataAER

! ****************************************************************
!                            
! CloseSaveDataAER - Open the data file and write header
!   Parameters :                                                  
!
! ****************************************************************
  
  SUBROUTINE CloseSaveDataAER ()
    ! modules
    USE {ROOT}_GlobalAER
    ! body
    CLOSE(fileidAER)
    CLOSE(86)       ! INTEGRATOR DIAGNOSTICS
    !
    DEALLOCATE(CSTAR)
    DEALLOCATE(GAMMA)
    DEALLOCATE(ORGANIC_SELECTION_INDICES)
    DEALLOCATE(RWORK)
    DEALLOCATE(IWORK)      
    DEALLOCATE(Y)

    END SUBROUTINE CloseSaveDataAER

! -----

end module {ROOT}_UtilAER
