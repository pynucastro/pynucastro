module network
  use physical_constants, only: ERG_PER_MeV
  
  implicit none

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 8

  ! Number of entries in reactvec returned by rate_evaluate
  integer, parameter :: nreactvec = 6
  
  ! Binding Energies Per Nucleon (MeV)
  double precision, dimension(nspec) :: ebind_per_nucleon

  ! Nucleon mass number A
  double precision, dimension(nspec) :: aion

  ! Nucleon atomic number Z
  double precision, dimension(nspec) :: zion
  
  ! Binding Energies (ergs)
  double precision, dimension(nspec) :: ebind

  ! Nuclides
  integer, parameter :: jn   = 1
  integer, parameter :: jp   = 2
  integer, parameter :: jhe4   = 3
  integer, parameter :: jc12   = 4
  integer, parameter :: jne20   = 5
  integer, parameter :: jne23   = 6
  integer, parameter :: jna23   = 7
  integer, parameter :: jmg23   = 8
  ! Energy Generation Rate
  integer, parameter :: jenuc   = 9

  ! Reactions
  integer, parameter :: k_c12_c12a_ne20   = 1
  integer, parameter :: k_c12_c12n_mg23   = 2
  integer, parameter :: k_c12_c12p_na23   = 3
  integer, parameter :: k_n_p             = 4
  integer, parameter :: k_na23_ne23       = 5
  integer, parameter :: k_ne23_na23       = 6

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

    ebind_per_nucleon(jn)   = 0.0d0
    ebind_per_nucleon(jp)   = 0.0d0
    ebind_per_nucleon(jhe4)   = 0.0d0
    ebind_per_nucleon(jc12)   = 0.0d0
    ebind_per_nucleon(jne20)   = 0.0d0
    ebind_per_nucleon(jne23)   = 0.0d0
    ebind_per_nucleon(jna23)   = 0.0d0
    ebind_per_nucleon(jmg23)   = 0.0d0
    
    aion(jn)   = 1.000000d+00
    aion(jp)   = 1.000000d+00
    aion(jhe4)   = 4.000000d+00
    aion(jc12)   = 1.200000d+01
    aion(jne20)   = 2.000000d+01
    aion(jne23)   = 2.300000d+01
    aion(jna23)   = 2.300000d+01
    aion(jmg23)   = 2.300000d+01

    zion(jn)   = 0.000000d+00
    zion(jp)   = 1.000000d+00
    zion(jhe4)   = 2.000000d+00
    zion(jc12)   = 6.000000d+00
    zion(jne20)   = 1.000000d+01
    zion(jne23)   = 1.000000d+01
    zion(jna23)   = 1.100000d+01
    zion(jmg23)   = 1.200000d+01

    do i = 1, nspec
      ebind(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do
  end subroutine init_network

end module network
