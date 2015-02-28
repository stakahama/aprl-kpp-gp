  module {ROOT}_Partition

  use {ROOT}_Precision, only: dp
  implicit none

contains

  ! ****************************************************************
  !
  ! partition - main partitioning program to be called from {ROOT}_Main.f90
  !   Parameters :
  !
  ! contains partitioning specific parameters, as well as the subroutine model_gas_aer
  !
  ! ****************************************************************

  subroutine partition(MCM_timestep_size, timesteps)
    ! This subroutine models gas-aerosol partitioning as described by Chen et al.
    ! (2010, Particle-Phase Chemistry of Secondary Organic Material: Modeled
    ! Compared to Measured O:C and H:C Elemental Ratios Provide Constraints)

    ! OUTPUT: the concentrations of different compounds will be saved in the arrays
    ! intermediate_CGAS and intermediate_CAER of size #compounds x #timesteps
    ! The last column of these arrays will be written into the global CGAS and CAER.

    ! This code calls a subroutine integrate() to compute one partitioning integration step
    ! CGAS_new and CAER_new will be computed by that subroutine. These variables are then, only for 
    ! the organic species, written into intermediate_CGAS and intermediate_CAER. At the end of the
    ! Partitioning step the global CGAS and CAER are updated and the code returns to main.

    ! modules
    use {ROOT}_Global, only: TEMP, NSPEC, CGAS => C, TIME, TEND, TSTART ! FB: to print the total TIME, TEND, and TSTART
    use {ROOT}_GlobalAER, only: CAER, VPAER, organic_selection_binary, d_ve_aerosols, nr_aerosol_particles, Diff_coeff, &
         displayed_not_enough_CGAS, displayed_not_enough_CAER, & ! FB: to count the error_messages
         f_a_Kn
    use {ROOT}_monitor, only: SPC_NAMES ! FB: to plot the spc_names

    ! local variables
    integer :: timesteps, current_step, j, i
    real (kind=dp), intent(in) :: MCM_timestep_size
    real (kind=dp), dimension(NSPEC, timesteps+1) :: intermediate_CGAS, intermediate_CAER       ! intermediate CAER and CGAS stock the integrated values in case multiple partitioning time steps are done within one MCM timestep
    real (kind=dp), dimension(NSPEC) :: CAER_new, CGAS_new, CAER_old, CGAS_old                  ! CAER_new, CAER_old, etc. are used for the starting and finishing values of the subroutine integrate() and are consequently written into intermediate_CAER and intermediate_CGAS
    !they might contain also the values for inorganic species. However, these are not propagated into intermediate_CAER

    real (kind=dp) :: PARTITION_timestep_size, rate, change

    ! body
    PARTITION_timestep_size = MCM_timestep_size/timesteps

    ! INTEGRATE
    intermediate_CAER = 0
    intermediate_CGAS = 0
    intermediate_CAER(:,1) = CAER
    intermediate_CGAS(:,1) = CGAS

    ! set error counter to zero for each MCM_timestep
    displayed_not_enough_CGAS = 0
    displayed_not_enough_CAER = 0

    do current_step = 1,timesteps
       ! subset the current, intermediate C and M to hand it over as CGAS_old and CAER_old to integrate to CGAS_new and CAER_new
       CAER_old = intermediate_CAER(:,current_step)
       CGAS_old = intermediate_CGAS(:,current_step)
       CAER_new = 0;
       CGAS_new = 0;

       call integrate(CGAS_old, CAER_old, CGAS_new, CAER_new, PARTITION_timestep_size, f_a_Kn, current_step)

       ! update only the organic compounds with the computed concentrations
       where(organic_selection_binary==1) intermediate_CAER(:,current_step+1) = CAER_new
       where(organic_selection_binary==1) intermediate_CGAS(:,current_step+1) = CGAS_new
       ! update the concentrations of the inorganic compounds with the preceeding value
       where(organic_selection_binary==0) intermediate_CAER(:,current_step+1) = CAER_old
       where(organic_selection_binary==0) intermediate_CGAS(:,current_step+1) = CGAS_old
    end do

    ! update concentrations in the global modules "Global" and "GlobalAER" to continue with MCM modelling
    CAER = intermediate_CAER(:,timesteps+1)
    CGAS = intermediate_CGAS(:,timesteps+1)
    !    do j = 1,NSPEC
    !      write(*,*) "current CGAS_old, planned CGAS:", CGAS(j), intermediate_CGAS(j,timesteps+1), 'for compound: ', SPC_NAMES(j)
    !      write(*,*) "current CGAS_old, planned CAER:", CAER(j), intermediate_CAER(j,timesteps+1), 'for compound: ', SPC_NAMES(j)
    !    end do

    ! write concentrations to file (Note: with the following loop we never write the initial concentration, this is happening in InitilaizeAER or in the timestep before)    do j = 1,timesteps
    do j = 1,timesteps
       write(86,987) TIME+(j)*PARTITION_timestep_size, (intermediate_CGAS(i,j+1), i=1,NSPEC)  ! implied loop to write all elements on an entire line of CGAS, this could also be done with two loops and write(*,*,advance="no")
       write(87,987) TIME+(j)*PARTITION_timestep_size, (intermediate_CAER(i,j+1), i=1,NSPEC)
    end do
987 FORMAT(1000ES19.11E3)  ! max 1000 species, then linebreak, can easily be increased here and in 986 (or use '<n>', but not supported by compiler gfortran)

  end subroutine partition





  subroutine integrate(CGAS_old, CAER_old, CGAS_new, CAER_new, dt, f_a_Kn, current_step)
    ! Note that integrate() computes (except for dlsode case), a change vector that includes inorganic species as well, these will just not be used by the calling routine

    ! modules
    use {ROOT}_Global, only: NSPEC, TIME, TEND, TSTART ! FB: to print the total TIME, TEND, and TSTART
    use {ROOT}_GlobalAER, only: VPAER, NorganicSPEC, organic_selection_binary, d_ve_aerosols, nr_aerosol_particles, Diff_coeff, &
         pi_constant, &
         molecular_masses, Avogadro, & ! !FB: to convert from molec/cm3 to µg/m
         displayed_not_enough_CGAS, displayed_not_enough_CAER, & ! FB: to count the error_messages
         organic_selection_indices ! FB: only used for dlsode (creation of special vector y (only containing organics, like gamma))
    use {ROOT}_monitor, only: SPC_NAMES ! FB: to plot the spc_names

    ! local variables
    real (kind=dp), intent(in) :: dt, f_a_Kn
    real (kind=dp), dimension(NSPEC) :: CGAS_new, CAER_new, CGAS_old, CAER_old
    integer, intent(in):: current_step
    character(len=100)::ERROR_MSG ! ERROR message to print

    real (kind=dp), dimension(NSPEC) :: factor, change
    real (kind=dp) :: sum_organic_CAER_old, sum_organic_CAER_old_microg_m3, mean_molecular_mass_of_organics
    integer :: ix

    ! LSODE PARAMETERS
    !external f1, jac1  ! NOTE (23.12.14): the declaration of f1 and jac1 seems not to be needed anymore in FORTRAN 90? (contradicting the DLSODE manual (lines 57 and 66 in opkdmain.f))
    integer i, iopt, istate, itask, itol, liw, lrw, mf, neq, ml, mu
    double precision atol, rtol, t, tout
    double precision, allocatable, dimension(:) :: rwork, y
    integer, allocatable, dimension(:) :: iwork
    real (kind=dp) :: dCGAS_i, dCAER_i ! needed for the check whether d CGAS/dt = -d CAER/dt for each compound

    ! A) INITIAL CHECK for negative concentrations compound and correct
    do ix = 1,NSPEC ! TODO: should I only do this check for organic compounds, or for all? only for organic: do index = 1,NorganicSPEC; ix = organic_selection_indices(index); ... leave the rest of the loop the same
       if (CAER_old(ix) < 0) then
          ! document by printing:
          ERROR_MSG = 'INTEGRATOR_ERROR: INITIAL_CHECK: correct negative CAER concentration: for compound: '
401       FORMAT(F6.1,'%. T=',ES9.3, ' ', A, I4, ' ', A, ' ', A, ES18.10E3, A, ES18.10E3, A, I3)
          write(*,401) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               ix, SPC_NAMES(ix), &
               ' CAER: ', CAER_old(ix), ' VPAER: ', VPAER(ix), '     at substep ', current_step
          write(88,401) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               ix, SPC_NAMES(ix), &
               ' CAER: ', CAER_old(ix), ' VPAER: ', VPAER(ix), '     at substep ', current_step
          ! and correct
          CAER_old(ix) = 0
       end if
       if (CGAS_old(ix) < 0) then
          ! document by printing:
          ERROR_MSG = 'INTEGRATOR_ERROR: INITIAL_CHECK: correct negative CGAS concentration: for compound: '
          write(*,401) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               ix, SPC_NAMES(ix), &
               ' CGAS: ', CGAS_old(ix), ' VPAER: ', VPAER(ix), '     at substep ', current_step
          write(88,401) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               ix, SPC_NAMES(ix), &
               ' CGAS: ', CGAS_old(ix), ' VPAER: ', VPAER(ix), '     at substep ', current_step
          ! and correct
          CGAS_old(ix) = 0
       end if
    end do      ! NOTE: without printing warnings the entire loop above could be written: "where(CAER_old < 0) CAER_old = 0" and "where(CGAS_old < 0) CGAS_old = 0"

    ! B) INTEGRATION
    change = 0
    ! calculate change in compounds
    sum_organic_CAER_old = 0   ! TODO: find a neater way, using only "organic_selection_binary & WHERE statement" instead of "organic_selection_indices & LOOP", NOTE ALREADY FOUND a way: "organic_selection_binary & LOOP"
    do ix = 1, NSPEC
       if (organic_selection_binary(ix)==1) then
          sum_organic_CAER_old = sum_organic_CAER_old + CAER_old(ix)
       end if
    end do
    mean_molecular_mass_of_organics = sum(CAER_old*molecular_masses)/sum(CAER_old)
    sum_organic_CAER_old_microg_m3 =sum_organic_CAER_old*(1000000_dp*mean_molecular_mass_of_organics/Avogadro)
    write(*,*) 'sum organic CAER before integrating [molecules/cm3] :', sum_organic_CAER_old, &
         'corresponds approx. to [microg/m3]: ', sum_organic_CAER_old_microg_m3

    factor = dt*2*pi_constant*nr_aerosol_particles*d_ve_aerosols*f_a_Kn*Diff_coeff;

    ! -------------------- DLSODE --------------------
    neq = NorganicSPEC*2
    mf = 21   ! defines the use of a dense Jacobian

    ! create state vector y = [CGAS; CAER]
    allocate(y(neq))
    do i=1,NorganicSPEC
       y(i) = CGAS_old(organic_selection_indices(i))
       y(i + NorganicSPEC) = CAER_old(organic_selection_indices(i))
    end do
    t = 0;
    tout = dt;
    ! other integrator parameters
    itol = 1 ! 1 or 2 according as ATOL (below) is a scalar or an array.
    !rtol = 0.0d0 ! Relative tolerance parameter (scalar).
    !atol = 1.0d+6 ! Absolute tolerance parameter (scalar or array).
    rtol = 1.0d-4 ! Relative tolerance parameter (scalar).
    atol = 1.0d-3 ! Absolute tolerance parameter (scalar or array).! The estimated local error in Y(i) will be controlled so as to be roughly less (in magnitude) than EWT(i) = RTOL*ABS(Y(i)) + ATOL
    itask = 1 ! Flag indicating the task DLSODE is to perform. Use ITASK = 1 for normal computation of output values of y at t = TOUT.
    istate = 1 ! Index used for input and output to specify the state of the calculation.
    iopt = 0 ! Flag indicating whether optional inputs are used (0: no, 1: yes)
    lrw = 22 + 9*neq + neq*neq ! Declared length of RWORK
    liw = 20 + neq             ! Declared length of IWORK
    allocate(rwork(lrw)) ! Real work array of length at least:    22 +  9*NEQ + NEQ**2     for MF = 21 or 22
    allocate(iwork(liw)) ! Integer work array of length at least: 20 + NEQ                 for MF = 21, 22, 24, or 25.        

    call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac1,mf)

    ! read out state vector y = [CGAS; CAER]
    do i=1,NorganicSPEC
       CGAS_new(organic_selection_indices(i)) = y(i)
       CAER_new(organic_selection_indices(i)) = y(i + NorganicSPEC)
    end do
    deallocate(y)
    deallocate(rwork)
    deallocate(iwork)

    ! check whether d CGAS/dt = -d CAER/dt 
    do i=1,NorganicSPEC
       dCGAS_i = (CGAS_old(organic_selection_indices(i)) - CGAS_new(organic_selection_indices(i)))
       dCAER_i = (CAER_old(organic_selection_indices(i)) - CAER_new(organic_selection_indices(i)))
       if (abs(dCGAS_i + dCAER_i)>atol) then
          ERROR_MSG = 'INTEGRATOR_ERROR: DLSODE OUTPUT: d CGAS/dt = -d CAER/dt is not fulfilled for compound: '
402       FORMAT(F6.1,'%. T=',ES9.3, ' ', A, I4, ' ', A, ' ', A, ES18.10E3, A, ES18.10E3, A, I3)
          write(*,402) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               organic_selection_indices(i), SPC_NAMES(organic_selection_indices(i)), &
               ' dCGAS/dt: ', dCGAS_i, ' dCAER/dt: ', dCAER_i, '     at substep ', current_step
          write(88,402) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
               organic_selection_indices(i), SPC_NAMES(organic_selection_indices(i)), &
               ' dCGAS/dt: ', dCGAS_i, ' dCAER/dt: ', dCAER_i, '     at substep ', current_step
       end if
    end do

    ! and carry on concentration of inorganic compounds
    do i=1,NSPEC
       if (organic_selection_binary(i) == 0) then
          CGAS_new(i) = CGAS_old(i)
          CAER_new(i) = CAER_old(i) ! or = 0
       end if
    end do

    ! to apply the check for negative concentrations below
    !change = CGAS_old-CGAS_new;    ! decide which one to use, they sould be equal after having passed the test above (check whether d CGAS/dt = -d CAER/dt)
    change = CAER_new - CAER_old;
    ! -------------------- END DLSODE --------------------

    ! -------------------- BEGIN CHECK --------------------
    ! C) FINAL CHECK for available amount of each compound and correct (only for forward and backward)
    
    ! check all the species
    do ix = 1,NSPEC ! TODO should I only do this check for organic compounds, or for all? (to change to only for organic: do index = 1,NorganicSPEC; ix = organic_selection_indices(index); ... leave the rest of the loop the same)
       ! check CAER
       if (CAER_old(ix) + change(ix) < 0) then
          ! display error message if it concerns organic species
          if (organic_selection_binary(ix) == 1 .AND. displayed_not_enough_CAER(ix) < 2) then
             ERROR_MSG = 'INTEGRATOR_ERROR: FINAL_CHECK: not enough CAER: for compound: '
403          FORMAT(F6.1,'%. T=',ES9.3, ' ', A, I4, ' ', A, ' ', A, ES18.10E3, A, ES18.10E3, A, I3)
             write(*,403) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
                  ix, SPC_NAMES(ix), &
                  'CAER_before: ', CAER_old(ix), ' intended change: ', change(ix), '     at substep ', current_step
             write(88,403) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
                  ix, SPC_NAMES(ix), &
                  'CAER_before: ', CAER_old(ix), ' intended change: ', change(ix), '     at substep ', current_step
          end if
          ! correct to zero
          change(ix) = - CAER_old(ix)
       end if
       ! check CGAS
       if (CGAS_old(ix) - change(ix) < 0) then
          ! display error message if it concerns organic species
          if (organic_selection_binary(ix) == 1 .AND. displayed_not_enough_CGAS(ix) < 2) then
             ERROR_MSG = 'INTEGRATOR_ERROR: FINAL_CHECK: not enough CGAS: for compound: '
             write(*,403) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
                  ix, SPC_NAMES(ix), &
                  'CGAS_before: ', CGAS_old(ix), ' intended change: ', change(ix), '     at substep ', current_step
             write(88,403) (TIME-TSTART)/(TEND-TSTART)*100, TIME, TRIM(ERROR_MSG), &
                  ix, SPC_NAMES(ix), &
                  'CGAS_before: ', CGAS_old(ix), ' intended change: ', change(ix), '     at substep ', current_step
          end if
          ! correct to zero
          change(ix) = CGAS_old(ix)
       end if
    end do ! TODO: without printing warnings the entire loop above could be written: "where(CAER_old + change < 0) change = -CAER_old" and "where(CAER_old + change < 0) change = CGAS_old"

    ! update concentrations
    CAER_new = CAER_old + change;
    CGAS_new = CGAS_old - change;
    ! if final check needed to correct the amount, set it manually to zero (numerically more stable)
    do ix = 1,NSPEC ! TODO, again: should I only do this check for organic compounds, or for all? only for organic: do index = 1,NorganicSPEC; ix = organic_selection_indices(index); ... leave the rest of the loop the same
       if (CAER_old(ix) + change(ix) < 0) then 
          CAER_new(ix) = 0 
       end if
       if (CGAS_old(ix) - change(ix) < 0) then
          CGAS_new(ix) = 0
       end if
    end do
    ! -------------------- END CHECK --------------------
  end subroutine integrate


  subroutine f1 (neq, t, y, ydot)
    ! modules
    use {ROOT}_GlobalAER, only: Cstar, gamma

    ! local variables
    integer neq
    double precision t, y, ydot
    dimension y(neq), ydot(neq)

    integer :: i, shift
    double precision sum_Caer

    shift = neq/2   ! shift between C_i and M_i of the same compound in the vector y

    ! precompute sum of Caer
    sum_Caer = 0
    do i=1,neq/2
       sum_Caer = sum_Caer + y(i+shift)
    end do

    ! compute Cgas_dot
    do i = 1,neq/2
       ydot(i) = gamma(i) * ( y(i) - y(i+shift) / sum_Caer * Cstar(i) )
    end do

    ! copy into Caer_dot
    do i = 1,neq/2
       ydot(shift + i) = - ydot(i)
    end do
    return
  end subroutine f1


  subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)

    !     DOCUMENTATION (l.262):
    !       " which provides df/dy by loading PD as follows:
    !        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
    !          the partial derivative of f(i) with respect to y(j).  (Ignore
    !          the ML and MU arguments in this case.)"

    ! modules
    use {ROOT}_GlobalAER, only: Cstar, gamma

    ! local variables
    integer neq, ml, mu, nrowpd
    double precision t, y, pd
    dimension y(neq), pd(nrowpd,neq)

    integer :: i, j, shift
    double precision sum_Caer

    shift = neq/2   ! shift between C_i and M_i of the same compound in the vector y

    ! precompute sum of Caer
    sum_Caer = 0
    do i=1,neq/2
       sum_Caer = sum_Caer + y(i+shift)
    end do

    ! fill the four different quadrants/sections of the Jacobian
    do i = 1,neq/2
       do j = 1,neq/2

          ! in comments: C = Cgas, M = Caer

          ! dCi'/dCj      = df(i)/dy(j)
          if (i==j) then
             pd(i,j) = gamma(i)
          end if ! otherwise pd(i,j)=0

          ! dCi'/dMj      = df(i)/dy(j+shift)
          if (i==j) then
             pd(i,j+shift) = - gamma(i) * Cstar(i) * ( 1/sum_Caer - y(i + shift)/(sum_Caer*sum_Caer))
          else
             pd(i,j+shift) = - gamma(i) * Cstar(i) * (            - y(i + shift)/(sum_Caer*sum_Caer))
          end if


          ! dMi'/dCj      = df(i+shift)/dy(j)
          if (i==j) then
             pd(i+shift,j) = - gamma(i)
          end if ! otherwise pd(i+shift,j)=0

          ! dMi'/dMj      = df(i+shift)/dy(j+shift)
          pd(i+shift,j+shift) = - pd(i,j+shift)
       end do
    end do
    return
  end subroutine jac1

end module {ROOT}_Partition
