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
    use {ROOT}_GlobalAER, only : CAER, VPAER, molecular_masses, &
         NorganicSPEC, cAER0_total, &
         Cstar, gamma, organic_selection_indices, &
         ! DLSODE
         iopt, istate, itask, itol, liw, lrw, mf, neq, ml, mu, &
         atol, rtol, rwork, y, iwork
    use {ROOT}_SIMPOLGroups, only: substruct_count, ngroups
    use {ROOT}_monitor, only: SPC_NAMES ! FB: to write the spc_names into the output file
    use simpol_module, only: simpolvp, press2conc

    ! local variables
    logical                     :: existp
    integer                     :: i, ix    
    !
    real(kind=dp)               :: alpha, lambda, Kn
    real(kind = dp), parameter  :: pi = 4_dp * atan(1.0_dp)
    real(kind = dp), parameter  :: Avogadro = 6.02214129E+23_dp ! [molecules/mol]
    real(kind=dp)               :: nr_aerosol_particles, d_ve_aerosols !FB
    real(kind=dp)               :: Diff_coeff(NSPEC)                   !FB
    real(kind=dp)               :: f_a_Kn                              ! FB: func(alpha, Knudsen-number)
    !
    real(dp)                    :: molefrac
    real(dp), dimension(NSPEC)  :: a0
    real(dp)                    :: CAER_total_in_g_per_m3, CAER_total_in_molecules_per_cm3, mean_molecular_mass_of_organics !FB: to convert initial amount of aerosol from µg/m3 to molecules/cm3
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
          vp_atm = simpolvp(TEMP, abundance)
          VPAER(i) = press2conc(TEMP, vp_atm) ! molecules/cm^3
       else
          VPAER(i) = 10.d12 ! no match; arbitrarily high number!-999.d0
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
    lrw = 22 + 9*neq + neq*neq ! Declared length of RWORK
    liw = 20 + neq             ! Declared length of IWORK
    mf = 21                    ! defines the use of a dense Jacobian
    allocate(rwork(lrw))       ! Real work array of length at least:    22 +  9*NEQ + NEQ**2     for MF = 21 or 22
    allocate(iwork(liw))       ! Integer work array of length at least: 20 + NEQ                 for MF = 21, 22, 24, or 25.
    allocate(y(neq))
    ! ------------------------------------------------------------    

    ! -------------------- INITIAL MOLE FRACTIONS --------------------
    a0 = 0.d0
    inquire(file="molefrac_init.txt", exist=existp)
    if (existp) then 
       open(55,file="molefrac_init.txt", status="old")
       read(55,*) !header/comment line
       do i = 1,NorganicSPEC
          read(55,*) ix, molefrac
          a0(ix) = molefrac
       end do
       close(55)
    else
       a0 = 1.d0
    endif
    ! ------------------------------------------------------------

    ! ------------------------------------------------------------    
    ! calculate initial CAER defined by the amount of CAER_total_in_g_per_m3     
    CAER_total_in_g_per_m3 =  cAER0_total ! [g/m3] defined by input text file (read in {ROOT}_Initialize.f90)
    ! old definition:       CAER_total_in_molecules_per_cm3 = NorganicSPEC*1E10 ! [molecules/cm3]     
    ! old definition corresponded to [g/m3]: 1598E-6 = 0.00159825 = 294*10^10*10^6*327.37/(6.022*10^23)     
    !
    ! convert to molecules/cm3
    mean_molecular_mass_of_organics = sum(a0*molecular_masses) ! [g/mol] to convert amount of initial total aerosol correctly from [g/m3] to [molecules/m3]
    write(*,*) "a0 weighted mean of molecular mass to transform I.C from microgram/m3 to molecules/cm3: ", &
         mean_molecular_mass_of_organics
    CAER_total_in_molecules_per_cm3 = CAER_total_in_g_per_m3/(1000000_dp*mean_molecular_mass_of_organics/Avogadro) ! [molecules/m3] = 100^3[cm3/m3][g/m3]/[g/molecule]
    CAER = CAER_total_in_molecules_per_cm3*a0    ! [molecules/cm3]
    ! ------------------------------------------------------------
    
    ! -------------------- prepare variables for integration (uses only the organic subset) --------------------
    allocate(Cstar(NorganicSPEC))
    allocate(gamma(NorganicSPEC))
    do i=1,NorganicSPEC
       Cstar(i) = VPAER(organic_selection_indices(i))
       gamma(i) = -2.d0*pi*nr_aerosol_particles*d_ve_aerosols*f_a_Kn*Diff_coeff(organic_selection_indices(i));
    end do
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
    use {ROOT}_GlobalAER, only : NorganicSPEC, organic_selection_indices, molecular_masses
    
    ! local variables
    integer :: ix
    real, dimension(:), allocatable    :: organics_mol_masses
    integer, dimension(:), allocatable :: organics

    include 'org_indices.f90' ! will define a vector 'organics' containing the indices of organics in the vectors "C0", "CAER", etc.
!!$      include 'org_molecular_masses.f90' ! will define a vector 'organics_mol_masses' containing the molecular masses of organics. It is ordered the same way as the indices in 'organics'.
    NorganicSPEC = size(organics)
    organic_selection_indices = organics

    molecular_masses = 0 ! inorganics will keep a molecular mass of 0

    do ix = 1,size(organics)
       ! save the molecular mass
       molecular_masses(organics(ix)) = organics_mol_masses(ix)
       !write(*,*) "read mol masses: ", molecular_masses(organics(ix)), "for comopound", organics(ix), ix, "th, organic compound"
    end do

  end subroutine read_in_subset_of_organics

end module {ROOT}_InitializeAER
