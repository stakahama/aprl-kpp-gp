  module {ROOT}_GlobalAER

  use {ROOT}_Precision, only: dp
  use {ROOT}_Parameters, only: NSPEC
  public
  save

  real(kind=dp) :: CAER(NSPEC)
  real(kind=dp) :: VPAER(NSPEC)

  real(kind=dp), parameter :: epsilon = .000001

  ! variables initialized by {ROOT}_Initialize.f90
  integer             :: partitioning_mode ! 2
  integer             :: absorptive_mode   ! 2
  real(kind=dp)       :: CAER0_total_microg_m3 = 0.d0  
  real(kind=dp)       :: CAER_total_microg_m3 = 0.d0        
  real(kind=dp)       :: CAER_total_molec_cm3 = 0.d0        
  real(kind=dp)       :: CAER_ghost_molec_cm3 = 0.d0
  integer             :: integratorcheck = 0
  logical             :: absorptivep=.TRUE.

  ! Define special variables to keep dlsode easy and organised (they contain only the organic species, unlike the fixed parameters above)
  real(kind=dp), dimension(:), allocatable :: Cstar  ! FB: VPAER but only containing organic species, will be used by DLSODE f-function
  real(kind=dp), dimension(:), allocatable :: gammaf ! FB: factor for partitioning of organic species, will be used by DLSODE f-function

  ! read in with read_in_subset_of_organics(), defined by org_indices.f90 and org_molecular_masses.f90
  integer                                  :: NorganicSPEC               !FB ! TODO: define whether this is needed after all (not if I'm only looping over the binary vector, and by defining CAER_total)
  integer, dimension(:), allocatable       :: organic_selection_indices  !FB: defines indices of organic compounds in "C0", "M0", etc.
  real(kind=dp), dimension(:), allocatable :: organic_molecular_masses
  real(kind=dp), dimension(:), allocatable :: xorgaer_init
  real(kind=dp), dimension(NSPEC)          :: orgmask
  real(kind=dp), dimension(NSPEC)          :: molecular_masses           !FB: molecular masses of organic compounds, (inorganic will have 0), to convert initial amount of aerosl from µg/m3 to molecules/cm3

  ! constants
  real(kind = dp), parameter  :: pi_constant = 4_dp * atan(1.0_dp)
  real(kind = dp), parameter  :: Avogadro = 6.02214129E+23_dp ! [molecules/mol]

  ! -------------------- DLSODE PARAMETERS --------------------    
  integer iopt, istate, itask, itol, liw, lrw, mf, neq, ml, mu
  double precision atol, rtol
  double precision, allocatable, dimension(:) :: rwork, y
  integer, allocatable, dimension(:)          :: iwork
  ! ------------------------------------------------------------

  real(kind=dp)                      :: minconc

  
  ! dynamic parameters
  real(kind=dp), parameter        :: alpha=1.d0            ! [ ] mass accommodation coefficient
  real(kind=dp), parameter        :: lambda=2.d-8          ! [m] mean free path
  real(kind=dp), parameter        :: rho_p=1.5d0           ! [g cm^-3] SOA density, Kostenidou et al. (2007)  
  real(kind=dp), dimension(NSPEC), parameter :: Diff_coeff=5.d-6      ! [m^2/s] diffusion coeff. of organic species in air
!!$  real(kind=dp)                   :: d_ve_aerosols ! [m] volume-equivalent diameter of the particles
  real(kind=dp)                   :: d_ve_seed     ! [m] volume-equivalent diameter of seed
  real(kind=dp)                   :: nr_aerosol_particles
  real(kind=dp)                   :: mean_molecular_mass_of_organics
  
!!$  nr_aerosol_particles, d_ve_seed are set in 'input_partitioning.txt'
  
!!$ Note mappings among arrays:  
!!$  x: CGAS(NSPEC), CAER(NSPEC), VPAER(NSPEC)
!!$  y: Cstar(NorganicSPEC), gammaf(NorganicSPEC)
!!$  do i=1,NorganicSPEC
!!$     iOrg=organic_selection_indices(i)
!!$     y(i)=fn(x(iOrg))
!!$  end do

contains

  subroutine print_organic_aerosol_mass()
  
    real(kind=dp) :: sum_Caer, org_mass

    sum_Caer = sum(CAER*orgmask)
    org_mass = calc_organic_aerosol_mass(CAER)
    write(*,*) 'sum organic CAER  [molecules/cm3] :', sum_Caer, &
         'corresponds approx. to [microg/m3]: ', org_mass

  end subroutine print_organic_aerosol_mass

  function calc_aerosol_conc(mass, mean_molec_mass) result(conc)
    ! mass is in microg/m^3
    ! conc is in molec/cm^3
    real(kind=dp), intent(in) :: mass
    real(kind=dp), intent(in) :: mean_molec_mass
    real(kind=dp)            :: conc

    ! [molecules/cm3] = [microg/m3][100^-6,m3/cm3*10^-6,g/microg][g/mole]^-1[molec/mole]
    conc = mass*1.d-12/mean_molec_mass*Avogadro

  end function calc_aerosol_conc

  function calc_organic_aerosol_mass(conc) result(org_mass)
    ! conc is in molec/cm^3
    ! org_mass is in microg/m^3
    real(kind=dp), dimension(:), intent(in) :: conc
    real(kind=dp)                           :: org_mass
    !
    if (size(conc) .eq. NSPEC) then
       org_mass = sum(conc*molecular_masses*orgmask)/Avogadro*1.d12
    else if (size(conc) .eq. NorganicSPEC) then
       org_mass = sum(conc*organic_molecular_masses)/Avogadro*1.d12   
    end if
    !
  end function calc_organic_aerosol_mass

  subroutine calc_aerosol_molefrac(caer,absorptivep,xorgaer)
    ! caer = molec/cm^3
    real(kind=dp), dimension(:), intent(in)  :: caer
    logical, intent(in)                      :: absorptivep
    real(kind=dp), dimension(:), intent(out) :: xorgaer
    ! local
    integer :: ix
    !
    if (absorptivep) then
       caer_total_molec_cm3 = sum(caer*orgmask)
       xorgaer = (/ (caer(organic_selection_indices(ix))/caer_total_molec_cm3, ix=1,NorganicSPEC) /)
    else 
       xorgaer = xorgaer_init
    end if
    !
  end subroutine calc_aerosol_molefrac

  function dp_from_mass(org_mass) result(diam)
    real(kind=dp), intent(in) :: org_mass
    real(kind=dp) :: diam
    ! org_mass = nr_aerosol_particles * rho_p * pi/6 * (diam**3 - diam_seed**3)
    diam = (org_mass / (nr_aerosol_particles*rho_p*pi_constant/6.d0) * 1.d-12 &
         + d_ve_seed**3.d0)**(1.d0/3.d0)
    end function dp_from_mass
  
  function calc_prefix(d_ve_aerosols,diff_coeff) result(gamma)
    real(kind=dp), intent(in) :: d_ve_aerosols
    real(kind=dp), intent(in) :: diff_coeff
    real(kind=dp)             :: gamma
    ! local vars:
    real(kind=dp)             :: Kn
    real(kind=dp)             :: f_a_Kn !FB: func(alpha, Knudsen-number)
    !
    Kn = 2.d0*lambda/d_ve_aerosols
    f_a_Kn = 0.75d0*alpha*(1 + Kn)/(Kn*Kn + Kn + 0.283d0*Kn*alpha + 0.75d0*alpha)
    !
    gamma = -2.d0*pi_constant*nr_aerosol_particles*d_ve_aerosols*f_a_Kn*diff_coeff
    !
  end function calc_prefix

end module {ROOT}_GlobalAER
