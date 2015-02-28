  module {ROOT}_GlobalAER

  use {ROOT}_Parameters, only: dp, NSPEC
  public
  save

  real (kind = dp), parameter :: Avogadro = 6.02214129E+23_dp ! [molecules/mol]
  real (kind = dp), parameter :: pi_constant = 4 * atan(1.0_dp)

  real(kind=dp) :: CAER(NSPEC)
  real(kind=dp) :: VPAER(NSPEC)

  ! variables initialized by {ROOT}_Initialize.f90
  integer             :: partition_substeps !FB
  real(kind=dp)       :: cAER0_total        !FB: initial total amount of aerosol present (default: 313E11 (NSPEC*E11), can be overwritten by input textfile)

  ! variables initialized by {ROOT}_InitializeAER.f90
  ! FIXED PARAMETERS (contain all the species, including anorganics)
  real(kind=dp) :: nr_aerosol_particles, d_ve_aerosols !FB
  real(kind=dp) :: Diff_coeff(NSPEC)                   !FB
  real(kind=dp) :: f_a_Kn                              ! FB: func(alpha, Knudsen-number)

  ! Define special variables to keep dlsode easy and organised (they contain only the organic species, unlike the fixed parameters above)
  real(kind=dp), dimension(:), allocatable :: Cstar ! FB: VPAER but only containing organic species, will be used by DLSODE f-function
  real(kind=dp), dimension(:), allocatable :: gamma ! FB: factor for partitioning of organic species, will be used by DLSODE f-function

  ! read in with read_in_subset_of_organics(), defined by org_indices.f90 and org_molecular_masses.f90
  integer             :: NorganicSPEC                              !FB ! TODO: define whether this is needed after all (not if I'm only looping over the binary vector, and by defining CAER_total)
  integer, dimension(:), allocatable :: organic_selection_indices  !FB: defines indices of organic compounds in "C0", "M0", etc.
  integer             :: organic_selection_binary(NSPEC)           !FB: defines (binary: 0 or 1) whether the compounds in "C0", "M0", etc. are organic and need to be modelled into the aerosol phase
  real(kind=dp)       :: molecular_masses(NSPEC)                   !FB: molecular masses of organic compounds, (inorganic will have 0), to convert initial amount of aerosl from µg/m3 to molecules/cm3

  ! variables used by {ROOT}_PartitionAER.f90
  integer :: displayed_not_enough_CAER(NSPEC)   !FB: flag to count how many times the error message on not enough CAER has already been shown for the given MCM_timestep
  integer :: displayed_not_enough_CGAS(NSPEC)   !FB: flag to count how many times the error message on not enough CAER has already been shown for the given MCM_timestep


end module {ROOT}_GlobalAER
