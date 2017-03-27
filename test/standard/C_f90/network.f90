module network
  use physical_constants, only: ERG_PER_MeV, me, mn, mp, CLIGHT, mev2gr
  
  implicit none

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 8

  integer, parameter :: nrat_tabular = 0

  integer, parameter :: nrat_reaclib = 5

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
  integer, parameter :: jn   = 1
  integer, parameter :: jp   = 2
  integer, parameter :: jhe4   = 3
  integer, parameter :: jc12   = 4
  integer, parameter :: jo16   = 5
  integer, parameter :: jne20   = 6
  integer, parameter :: jna23   = 7
  integer, parameter :: jmg23   = 8

  integer :: net_ienuc = nspec + 1

  ! Reactions
  integer, parameter :: k_c12_c12a_ne20   = 1
  integer, parameter :: k_c12_c12n_mg23   = 2
  integer, parameter :: k_c12_c12p_na23   = 3
  integer, parameter :: k_c12_ag_o16   = 4
  integer, parameter :: k_n_p   = 5

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

    ebind_per_nucleon(jn)   = 0.00000000000000d+00
    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jna23)   = 8.11149300000000d+00
    ebind_per_nucleon(jmg23)   = 7.90110400000000d+00
    
    aion(jn)   = 1.00000000000000d+00
    aion(jp)   = 1.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jna23)   = 2.30000000000000d+01
    aion(jmg23)   = 2.30000000000000d+01

    zion(jn)   = 0.00000000000000d+00
    zion(jp)   = 1.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jne20)   = 1.00000000000000d+01
    zion(jna23)   = 1.10000000000000d+01
    zion(jmg23)   = 1.20000000000000d+01

    nion(jn)   = 1.00000000000000d+00
    nion(jp)   = 0.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jna23)   = 1.20000000000000d+01
    nion(jmg23)   = 1.10000000000000d+01
    
    do i = 1, nspec
       ebind(i) = ebind_per_nucleon(i) * aion(i) 
    end do

    ! Set the mass energy in erg
    mion(:) = (nion(:) * mn + zion(:) * (mp + me) &
         - ebind(:) * mev2gr)*CLIGHT*CLIGHT

  end subroutine init_network

end module network
