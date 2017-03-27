module network
  use physical_constants, only: ERG_PER_MeV, me, mn, mp, CLIGHT, mev2gr
  
  implicit none

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 2

  integer, parameter :: nrat_tabular = 0

  integer, parameter :: nrat_reaclib = 2

  ! Number of entries in reactvec returned by rate_evaluate
  integer, parameter :: nreactvec = 6

  ! Binding Energies Per Nucleon (MeV)
  double precision, dimension(nspec) :: ebind_per_nucleon

  ! Nucleon mass number A
  double precision, dimension(nspec) :: aion

  ! Nucleon atomic number Z
  double precision, dimension(nspec) :: zion

  ! Nucleon atomic number N
  double precision, dimension(nspec) :: nion

  ! Nucleon mass energy
  double precision, dimension(nspec) :: mion
  
  ! Binding Energies (ergs)
  double precision, dimension(nspec) :: ebind

  ! Nuclides
  integer, parameter :: jhe4   = 1
  integer, parameter :: jc12   = 2

  integer :: net_ienuc = nspec + 1

  ! Reactions
  integer, parameter :: k_c12_gaa_he4   = 1
  integer, parameter :: k_he4_aag_c12   = 2

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

contains

  subroutine init_network()
    integer :: i

    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01

    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00

    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    
    do i = 1, nspec
       ebind(i) = ebind_per_nucleon(i) * aion(i) 
    end do

    ! Set the mass energy in erg
    mion(:) = (nion(:) * mn + zion(:) * (mp + me) &
         - ebind(:) * mev2gr)*CLIGHT*CLIGHT

  end subroutine init_network

end module network
