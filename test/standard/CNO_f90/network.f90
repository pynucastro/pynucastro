module network
  use physical_constants, only: ERG_PER_MeV, me, mn, mp, CLIGHT, mev2gr
  
  implicit none

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 9

  integer, parameter :: nrat_tabular = 0

  integer, parameter :: nrat_reaclib = 8

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
  integer, parameter :: jp   = 1
  integer, parameter :: jhe4   = 2
  integer, parameter :: jc12   = 3
  integer, parameter :: jc13   = 4
  integer, parameter :: jn13   = 5
  integer, parameter :: jn14   = 6
  integer, parameter :: jn15   = 7
  integer, parameter :: jo14   = 8
  integer, parameter :: jo15   = 9

  integer :: net_ienuc = nspec + 1

  ! Reactions
  integer, parameter :: k_c12_pg_n13   = 1
  integer, parameter :: k_n13_c13   = 2
  integer, parameter :: k_c13_pg_n14   = 3
  integer, parameter :: k_n14_pg_o15   = 4
  integer, parameter :: k_o15_n15   = 5
  integer, parameter :: k_n15_pa_c12   = 6
  integer, parameter :: k_n13_pg_o14   = 7
  integer, parameter :: k_o14_n14   = 8

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

    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jc13)   = 7.46984900000000d+00
    ebind_per_nucleon(jn13)   = 7.23886300000000d+00
    ebind_per_nucleon(jn14)   = 7.47561400000000d+00
    ebind_per_nucleon(jn15)   = 7.69946000000000d+00
    ebind_per_nucleon(jo14)   = 7.05230100000000d+00
    ebind_per_nucleon(jo15)   = 7.46369200000000d+00
    
    aion(jp)   = 1.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jc13)   = 1.30000000000000d+01
    aion(jn13)   = 1.30000000000000d+01
    aion(jn14)   = 1.40000000000000d+01
    aion(jn15)   = 1.50000000000000d+01
    aion(jo14)   = 1.40000000000000d+01
    aion(jo15)   = 1.50000000000000d+01

    zion(jp)   = 1.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jc13)   = 6.00000000000000d+00
    zion(jn13)   = 7.00000000000000d+00
    zion(jn14)   = 7.00000000000000d+00
    zion(jn15)   = 7.00000000000000d+00
    zion(jo14)   = 8.00000000000000d+00
    zion(jo15)   = 8.00000000000000d+00

    nion(jp)   = 0.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jc13)   = 7.00000000000000d+00
    nion(jn13)   = 6.00000000000000d+00
    nion(jn14)   = 7.00000000000000d+00
    nion(jn15)   = 8.00000000000000d+00
    nion(jo14)   = 6.00000000000000d+00
    nion(jo15)   = 7.00000000000000d+00
    
    do i = 1, nspec
       ebind(i) = ebind_per_nucleon(i) * aion(i) 
    end do

    ! Set the mass energy in erg
    mion(:) = (nion(:) * mn + zion(:) * (mp + me) &
         - ebind(:) * mev2gr)*CLIGHT*CLIGHT

  end subroutine init_network

end module network
