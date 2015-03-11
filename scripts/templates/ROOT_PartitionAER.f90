  module {ROOT}_Partition

  use {ROOT}_Precision, only: dp
  implicit none
  private
  public partition

contains

  ! *************************************************************************** !  
  !
  ! partition - main partitioning program to be called from {ROOT}_Main.f90
  !   Parameters :
  !
  ! *************************************************************************** !  


  ! *************************************************************************** !
  subroutine partition(t, tout)

    ! modules
    use {ROOT}_Parameters, only: NSPEC
    use {ROOT}_Global, only: CGAS => C
    use {ROOT}_GlobalAER, only: CAER, print_organic_aerosol_mass

    ! bound variables
    real(kind=dp), intent(in) :: t, tout
    
    ! local variables
    real(kind=dp), dimension(NSPEC) :: CAER_new, CGAS_new
    
    call print_organic_aerosol_mass()

    CGAS_new = CGAS;
    CAER_new = CAER;

    call integrate(t, tout, CGAS, CAER, CGAS_new, CAER_new)

    CGAS = CGAS_new
    CAER = CAER_new
    
987 FORMAT(1000ES19.11E3)  ! max 1000 species, then linebreak, can easily be increased here and in 986 (or use '<n>', but not supported by compiler gfortran)

  end subroutine partition
  ! *************************************************************************** !

  ! *************************************************************************** !  
  subroutine integrate(TINIT, TNEXT, CGAS_old, CAER_old, CGAS_new, CAER_new)
    ! Only calculates change for organic species
    ! Partitioning does not affect concentrations of inorganic species

    ! -------------------- variables --------------------
    ! modules
    use {ROOT}_Parameters, only: NSPEC
    use {ROOT}_Global, only: TIME, TEND, TSTART ! FB: to print the total TIME, TEND, and TSTART
    use {ROOT}_GlobalAER, only: NorganicSPEC, &
         integratorcheck, Avogadro, &           
         organic_selection_indices, & 
         minconc, absorptivep, molecular_masses, &
         CAER_total_molec_cm3, CAER_total_microg_m3, CAER0_total_microg_m3, &
         ! DLSODE (values set in _Initialize_AER.f90)
         iopt, istate, itask, itol, liw, lrw, mf, neq, ml, mu, &
         atol, rtol, rwork, y, iwork
    use {ROOT}_monitor, only: SPC_NAMES ! FB: to plot the spc_names
    ! bound variables
    real(kind=dp), intent(in)                    :: TINIT, TNEXT
    real(kind=dp), dimension(NSPEC), intent(in)  :: CGAS_old, CAER_old
    real(kind=dp), dimension(NSPEC), intent(out) :: CGAS_new, CAER_new
    ! local variables
    real(kind=dp)     :: T, TOUT 
    character(len=100):: ERROR_MSG                   ! ERROR message to print
!!$    real (kind=dp) :: sum_organic_CAER_old, sum_organic_CAER_old_microg_m3, mean_molecular_mass_of_organics
    integer           :: i, iAER, iOrg
    real(kind=dp)     :: conversionf, change, dGAS, dAER ! needed for the check whether d CGAS/dt = -d CAER/dt for each compound
    ! LSODE PARAMETERS
    !external f1, jac1  ! NOTE (23.12.14): the declaration of f1 and jac1 seems not to be needed anymore in FORTRAN 90? (contradicting the DLSODE manual (lines 57 and 66 in opkdmain.f))
    ! ------------------------------------------------------------

    ! -------------------- integrate --------------------
    T = 0.0_dp
    TOUT = TNEXT-TINIT ! DT
    istate = 1
    y = 0.0_dp
    rwork = 0.0_dp
    iwork = 0
    ! create state vector y = [CGAS; CAER]
    CAER_total_molec_cm3 = 0 ! for integrator
    CAER_total_microg_m3 = 0 ! for absorptivep
    conversionf = 1.d12/Avogadro
    do i=1,NorganicSPEC
       iAER = i+NorganicSPEC
       iOrg = organic_selection_indices(i)
       y(i) = max(CGAS_old(iOrg),minconc)     ! this stabilizes the solution 
       y(iAER) = max(CAER_old(iOrg),minconc)  ! this stabilizes the solution 
       CAER_total_molec_cm3 = CAER_total_molec_cm3 + y(iAER)
       CAER_total_microg_m3 = CAER_total_microg_m3 + y(iAER)*molecular_masses(iOrg)*conversionf
    end do
    !
    if (.NOT. absorptivep .AND. CAER_total_microg_m3 .GE. CAER0_total_microg_m3) then
       absorptivep = .TRUE.
    end if
    !
    ! integrate
    if (mf .eq. 21) then
       call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac1,mf)
    else
       call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,dummy,mf)
    endif
    write(*,*) "Partitioning LSODE ISTATE=", istate
    ! ------------------------------------------------------------

    ! -------------------- postprocess --------------------
    ! read and assign state vector y = [CGAS; CAER]
    do i=1,NorganicSPEC
       iAER = i+NorganicSPEC
       iOrg = organic_selection_indices(i)
       if (integratorcheck > 0) then
          ! check whether d CGAS/dt = -d CAER/dt
          dGAS = y(i)-CGAS_old(iOrg)
          dAER = y(iAER)-CAER_old(iOrg)
          if (abs((dGAS+dAER)) > atol) then
             ERROR_MSG = 'INTEGRATOR_ERROR: DLSODE OUTPUT: d CGAS/dt = -d CAER/dt is not fulfilled for compound: '
             write(88,402) (T-TSTART)/(TEND-TSTART)*100, T, TRIM(ERROR_MSG), &
                  iOrg, SPC_NAMES(iOrg), &
                  ' dCGAS/dt: ', dGAS, ' dCAER/dt: ', dAER
          endif
       endif
       dGAS = max(y(i),minconc)-CGAS_old(iOrg)
       dAER = max(y(iAER),minconc)-CAER_old(iOrg)
       change = sign(min(abs(dGAS),abs(dAER)),dGAS)
       CGAS_new(iOrg) = CGAS_old(iOrg)+change
       CAER_new(iOrg) = CAER_old(iOrg)-change
    end do
    ! ------------------------------------------------------------

402 FORMAT(F6.1,'%. T=',ES9.3, ' ', A, I4, ' ', A, ' ', A, ES18.10E3, A, ES18.10E3, A, I3)
  end subroutine integrate
  ! --------------------------------------------------------------------------- !

subroutine dummy (neq, t, y, ml, mu, pd, nrowpd)

    ! local variables
    integer neq, ml, mu, nrowpd
    double precision t, y, pd
    dimension y(neq), pd(nrowpd,neq)

    ! body: pass
end subroutine dummy
  

  ! -------------------- routines required by LSODE --------------------

  ! *************************************************************************** !
  subroutine f1 (neq, t, y, ydot)
    ! modules
    use {ROOT}_GlobalAER, only: Cstar, gammaf, sum_Caer => CAER_total_molec_cm3, absorptivep

    ! local variables
    integer          :: neq
    double precision :: t, y, ydot
    dimension y(neq), ydot(neq)

    integer          :: i, nspec
    double precision :: x(neq)

    nspec = neq/2   ! shift between C_i and M_i of the same compound in the vector y

    ! calculate mole fraction
    if (absorptivep) then
       do i = 1,nspec
          x(i) = y(i+nspec)/sum_Caer ! mole fraction
       end do
    else
       x = 1.d0
    end if

    ! compute Cgas_dot; Caer_dot
    do i = 1,nspec
       ydot(i) = gammaf(i) * ( y(i) - x(i) * Cstar(i) )
       ydot(i+nspec) = - ydot(i)    
    end do

    return
  end subroutine f1
  ! *************************************************************************** !

  ! *************************************************************************** !
  subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)

    !     DOCUMENTATION (l.262):
    !       " which provides df/dy by loading PD as follows:
    !        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
    !          the partial derivative of f(i) with respect to y(j).  (Ignore
    !          the ML and MU arguments in this case.)"

    ! modules
    use {ROOT}_GlobalAER, only: Cstar, gammaf, sum_Caer => CAER_total_molec_cm3

    ! local variables
    integer neq, ml, mu, nrowpd
    double precision t, y, pd
    dimension y(neq), pd(nrowpd,neq)

    integer :: i, j, nspec
!!$    double precision sum_Caer

    nspec = neq/2   ! shift between C_i and M_i of the same compound in the vector y

    ! fill the four different quadrants/sections of the Jacobian
    do i = 1,nspec
       do j = 1,nspec

          ! in comments: C = Cgas, M = Caer

          ! dCi'/dCj      = df(i)/dy(j)
          if (i==j) then
             pd(i,j) = gammaf(i)
          end if ! otherwise pd(i,j)=0

          ! dCi'/dMj      = df(i)/dy(j+nspec)
          if (i==j) then
             pd(i,j+nspec) = - gammaf(i) * Cstar(i) * ( 1/sum_Caer - y(i+nspec)/(sum_Caer*sum_Caer))
          else
             pd(i,j+nspec) = - gammaf(i) * Cstar(i) * (            - y(i+nspec)/(sum_Caer*sum_Caer))
          end if


          ! dMi'/dCj      = df(i+nspec)/dy(j)
          if (i==j) then
             pd(i+nspec,j) = - gammaf(i)
          end if ! otherwise pd(i+nspec,j)=0

          ! dMi'/dMj      = df(i+nspec)/dy(j+nspec)
          pd(i+nspec,j+nspec) = - pd(i,j+nspec)
       end do
    end do
    
    return
  end subroutine jac1
  
end module {ROOT}_Partition
