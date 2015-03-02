  module {ROOT}_GlobalAER

  use {ROOT}_Precision, only: dp
  use {ROOT}_Parameters, only: NSPEC
  public
  save

  real(kind=dp) :: CAER(NSPEC)
  real(kind=dp) :: VPAER(NSPEC)

  ! variables initialized by {ROOT}_Initialize.f90
  integer             :: partition_on !FB
  real(kind=dp)       :: cAER0_total        !FB: initial total amount of aerosol present (default: 313E11 (NSPEC*E11), can be overwritten by input textfile)
  integer             :: integratorcheck = 0

  ! Define special variables to keep dlsode easy and organised (they contain only the organic species, unlike the fixed parameters above)
  real(kind=dp), dimension(:), allocatable :: Cstar  ! FB: VPAER but only containing organic species, will be used by DLSODE f-function
  real(kind=dp), dimension(:), allocatable :: gammaf ! FB: factor for partitioning of organic species, will be used by DLSODE f-function

  ! read in with read_in_subset_of_organics(), defined by org_indices.f90 and org_molecular_masses.f90
  integer                            :: NorganicSPEC               !FB ! TODO: define whether this is needed after all (not if I'm only looping over the binary vector, and by defining CAER_total)
  integer, dimension(:), allocatable :: organic_selection_indices  !FB: defines indices of organic compounds in "C0", "M0", etc.
  real(kind=dp), dimension(NSPEC)    :: molecular_masses           !FB: molecular masses of organic compounds, (inorganic will have 0), to convert initial amount of aerosl from µg/m3 to molecules/cm3

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

contains

  subroutine print_organic_aerosol_mass()
  
    real(kind=dp) :: sum_aer, sum_prod, mean_molec_mass, sum_aer_microg_m3
    integer :: ix, iOrg

    sum_aer = 0
    sum_prod = 0
    do ix = 1, NorganicSPEC
       iOrg = organic_selection_indices(ix)
       sum_aer = sum_aer + CAER(iOrg)
       sum_prod = sum_prod + (CAER(iOrg)*molecular_masses(iOrg))
    end do
    mean_molec_mass = sum_prod/sum_aer
    sum_aer_microg_m3 =sum_aer*(1000000_dp*mean_molec_mass/Avogadro)
    write(*,*) 'sum organic CAER  [molecules/cm3] :', sum_aer, &
         'corresponds approx. to [microg/m3]: ', sum_aer_microg_m3

  end subroutine print_organic_aerosol_mass

end module {ROOT}_GlobalAER
