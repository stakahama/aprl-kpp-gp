  module {ROOT}_InitializeAER
  implicit none
  private
  public InitializeAER

contains

  subroutine InitializeAER()
    ! import modules
    use {ROOT}_Precision, only : dp
    use {ROOT}_Parameters, only : NSPEC
    use {ROOT}_Global, only : TEMP, CFACTOR, CGAS => C, TSTART
    use {ROOT}_GlobalAER, only : CAER, VPAER, NorganicSPEC, &
         molecular_masses, orgmask, organic_selection_indices, &
         CAER0_total_microg_m3, CAER_total_microg_m3, CAER_total_molec_cm3, &
         Cstar, gammaf, pi_constant, &
         calc_aerosol_conc, absorptivep, absorptive_mode, &
         xorgaer_init, epsilon, CAER_ghost_molec_cm3, Avogadro, &
         ! DLSODE
         iopt, istate, itask, itol, liw, lrw, mf, neq, ml, mu, &
         atol, rtol, rwork, y, iwork
    use {ROOT}_SIMPOLGroups, only: substruct_count, ngroups
    use {ROOT}_monitor, only: SPC_NAMES ! FB: to write the spc_names into the output file
    use simpol_module, only: simpolvp, press2conc

    ! local variables
    logical                     :: existp
    integer                     :: i, ix, iOrg
    !
    real(kind=dp)               :: alpha, lambda, Kn
    real(kind=dp)               :: nr_aerosol_particles, d_ve_aerosols !FB
    real(kind=dp)               :: Diff_coeff(NSPEC)                   !FB
    real(kind=dp)               :: f_a_Kn                              ! FB: func(alpha, Knudsen-number)
    !
    real(dp)                    :: molefrac
    real(dp), dimension(NSPEC)  :: a0
    real(dp)                    :: mean_molecular_mass_of_organics !FB: to convert initial amount of aerosol from µg/m3 to molecules/cm3
    !
    integer, dimension(ngroups) :: abundance
    real(dp)                    :: vp_atm=0.d0

    ! body
    ! -------------------- FIXED PARAMETERS --------------------
    Diff_coeff = 5.d-6          ! [m^2/s] diffusion coeff. of organic species in air
    nr_aerosol_particles = 5.d9 ! [m^-3] particle number concentration
    d_ve_aerosols = 100.d-9     ! [m] volume-equivalent diameter of the particles
    alpha = 1.d0
    lambda = 2.d-8              ! [m] mean free path
    Kn = 2.d0*lambda/d_ve_aerosols
    f_a_Kn = 0.75d0*alpha*(1 + Kn)/(Kn*Kn + Kn + 0.283d0*Kn*alpha + 0.75d0*alpha)
    ! ------------------------------------------------------------    

    ! -------------------- VPAER from SIMPOL calculations --------------------
    do i=1,NSPEC
       abundance = substruct_count(i)
       if (any(abundance(2:) .gt. 0.)) then
          VPAER(i) = simpolvp(TEMP, abundance) ! [atm]
       else
          VPAER(i) = 10.d0                     ! no match; arbitrarily high number
       end if
    end do
    ! ------------------------------------------------------------

    ! ------------------------------------------------------------
    ! CAER (zero for inorganics, as we model only organic partitioning)
    ! define subset of organics
    molecular_masses = 0
    call read_in_subset_of_organics()
    ! assigns values to:
    !    NorganicSPEC
    !    organic_selection_indices
    !    molecular_masses
    ! ------------------------------------------------------------    

    ! -------------------- DLSODE PARAMETERS --------------------    
    neq = NorganicSPEC*2
    itol = 1                   ! 1 or 2 according as ATOL (below) is a scalar or an array.
    !rtol = 0.0d0 ! Relative tolerance parameter (scalar).
    !atol = 1.0d+6 ! Absolute tolerance parameter (scalar or array).
    rtol = 1.0d-4              ! Relative tolerance parameter (scalar).
    atol = 1.0d-3              ! Absolute tolerance parameter (scalar or array).! The estimated local error in Y(i) will be controlled so as to be roughly less (in magnitude) than EWT(i) = RTOL*ABS(Y(i)) + ATOL
    itask = 1                  ! Flag indicating the task DLSODE is to perform. Use ITASK = 1 for normal computation of output values of y at t = TOUT.
    istate = 1                 ! Index used for input and output to specify the state of the calculation.
    iopt = 0                   ! Flag indicating whether optional inputs are used (0: no, 1: yes)
!!$    mf = 10
!!$    mf = 21
!!$    mf = 22
    !
    select case (mf)              ! Read in from {ROOT}_Initialize.f90
    case (10)                     ! Nonstiff
       lrw = 20+16*neq            ! Declared length of RWORK
       liw = 20                   ! Declared length of IWORK
       !
    case (21)                     ! User-supplied dense Jacobian
       lrw = 22 + 9*neq + neq*neq ! Declared length of RWORK
       liw = 20 + neq             ! Declared length of IWORK
       ! 
    case (22)                     ! Internally generated Jacobian
       lrw = 22 + 9*neq + neq*neq ! Declared length of RWORK
       liw = 20 + neq             ! Declared length of IWORK
    end select
    !
    allocate(rwork(lrw))       ! Real work array of length at least:    22 +  9*NEQ + NEQ**2     for MF = 21 or 22
    allocate(iwork(liw))       ! Integer work array of length at least: 20 + NEQ                 for MF = 21, 22, 24, or 25.
    allocate(y(neq))
    ! ------------------------------------------------------------    
    
    ! -------------------- prepare variables for integration (uses only the organic subset) --------------------
    allocate(Cstar(NorganicSPEC))
    allocate(gammaf(NorganicSPEC))
    do i=1,NorganicSPEC
       iOrg = organic_selection_indices(i)
       Cstar(i) = press2conc(TEMP,VPAER(iOrg)) ! [molec/cm^3]
       gammaf(i) = -2.d0*pi_constant*nr_aerosol_particles*d_ve_aerosols*f_a_Kn*Diff_coeff(iOrg);
    end do
    ! ------------------------------------------------------------    

    ! -------------------- INITIAL MOLE FRACTIONS --------------------
    allocate(xorgaer_init(NorganicSPEC))
    a0 = 1.d0
    xorgaer_init = 1.d0
    inquire(file="molefrac_init.txt", exist=existp)
    if (existp) then 
       open(55,file="molefrac_init.txt", status="old")
       read(55,*) !header/comment line (skip)
       do i = 1,NorganicSPEC
          read(55,*) ix, molefrac
          a0(ix) = molefrac
       end do
       close(55)
       xorgaer_init = a0
       !
       if(all(a0 .le. epsilon)) then
          mean_molecular_mass_of_organics = 1.d0          
!!$       else if(all(abs(a0-1.d0) .le. epsilon)) then
!!$          mean_molecular_mass_of_organics = 1.d0          
       else
          a0 = a0/sum(a0) ! ensure mole fractions == 1
          mean_molecular_mass_of_organics = sum(a0*molecular_masses*orgmask)
       end if
       !
       if (absorptive_mode .eq. 0) then
          absorptivep = .TRUE.
          write(*,*) "a0 weighted mean of molecular mass to transform I.C from microgram/m3 to molecules/cm3: ", &
               mean_molecular_mass_of_organics
          CAER_total_molec_cm3 = calc_aerosol_conc(CAER0_total_microg_m3, mean_molecular_mass_of_organics)
       else if (absorptive_mode .gt. 0) then
          absorptivep = .FALSE.
          CAER_total_molec_cm3 = 0.d0
       end if
       CAER = CAER_total_molec_cm3*a0    ! [molecules/cm3]
    else
       ! self-starting
       absorptivep = .TRUE.
       CAER = 0.d0
       CAER_total_molec_cm3 = 0.d0
       do i=1,NorganicSPEC
          iOrg = organic_selection_indices(i)
          CAER(iOrg) = CAER0_total_microg_m3*CGAS(iOrg)/CSTAR(i)/molecular_masses(iOrg)*Avogadro*1.d-12
          CAER_total_molec_cm3 = CAER_total_molec_cm3 + CAER(iOrg)
       end do
       mean_molecular_mass_of_organics = sum(CAER*molecular_masses*ORGMASK)/CAER_total_molec_cm3
       if(absorptive_mode .eq. 0) then
          CAER_ghost_molec_cm3 = (CAER0_total_microg_m3/mean_molecular_mass_of_organics*Avogadro*1.d-12)-CAER_total_molec_cm3
       else if(absorptive_mode .gt. 0) then
          CAER_ghost_molec_cm3 = 0.d0
       end if
       xorgaer_init = CAER/CAER_total_molec_cm3
       !
    endif
    !
    ! ------------------------------------------------------------

    ! -------------------- open files and overwrite them completely --------------------
    inquire(file="output_ERRORS.txt", exist=existp)
    if (existp) then
       open(88, file="output_ERRORS.txt", status="old", position="rewind", action="write")
    else
       open(88, file="output_ERRORS.txt", status="new", action="write")
    end if
    ! ------------------------------------------------------------
    
986 FORMAT(1000(1X,A))  ! max 1000 species, then linebreak, can easily be increased here and in 987 (or use '<n>', but not supported by compiler gfortran)
    ! 988 same as 987 in {ROOT}_PartitionAER.f90
988 FORMAT(1000ES19.11E3)  ! max 1000 species, then linebreak, can easily be increased here and in 986 (or use '<n>', but not supported by compiler gfortran)

  end subroutine InitializeAER

  subroutine read_in_subset_of_organics()
    ! A vector "organics" with the indices of these organic compounds has been created before compilation with a python script.
    ! "organic_selection_indices" is saved in a file org_indices.f90.

    ! OUTPUT: the vectors "organic_selection_indices","organic_selection_binary", and molecular_masses will be overwritten
    use {ROOT}_Parameters, only : NSPEC
    use {ROOT}_GlobalAER, only : NorganicSPEC, &
         organic_selection_indices, organic_molecular_masses, &
         molecular_masses, orgmask
    
    ! local variables
    integer :: ix
    integer, dimension(:), allocatable :: organics
    real, dimension(:), allocatable    :: organics_mol_masses

    include 'org_indices.f90' ! will define a vector 'organics' containing the indices of organics in the vectors "C0", "CAER", etc.

    !
    NorganicSPEC = size(organics)
    allocate(organic_selection_indices(NorganicSPEC))
    allocate(organic_molecular_masses(NorganicSPEC))
    organic_selection_indices = organics
    organic_molecular_masses = organics_mol_masses

    orgmask = 0.d0
    molecular_masses = 0.d0 ! inorganics will keep a molecular mass of 0
    do ix = 1,NorganicSPEC
       ! save the molecular mass
       molecular_masses(organics(ix)) = organics_mol_masses(ix)
       orgmask(organics(ix)) = 1.d0
       !write(*,*) "read mol masses: ", molecular_masses(organics(ix)), "for comopound", organics(ix), ix, "th, organic compound"
    end do
    
  end subroutine read_in_subset_of_organics

end module {ROOT}_InitializeAER
