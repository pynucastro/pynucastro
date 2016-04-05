module network
  use physical_constants, only: ERG_PER_MeV
  
  implicit none

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 7

  ! Binding Energies Per Nucleon (MeV)
  double precision, dimension(nspec) :: ebind_per_nucleon

  ! Nucleon mass number A
  double precision, dimension(nspec) :: aion

  ! Nucleon atomic number Z
  double precision, dimension(nspec) :: zion
  
  ! Binding Energies (ergs)
  double precision, dimension(nspec) :: ebind

  ! Nuclides
  integer :: jn   = 1
  integer :: jp   = 2
  integer :: jhe4   = 3
  integer :: jc12   = 4
  integer :: jne20   = 5
  integer :: jna23   = 6
  integer :: jmg23   = 7
  ! Energy Generation Rate
  integer :: jenuc   = 8

  ! Reactions
  integer :: k_c12_c12a_ne20   = 1
  integer :: k_c12_c12n_mg23   = 2
  integer :: k_c12_c12p_na23   = 3
  integer :: k_n_p   = 4
  
contains

  subroutine init_network()
    integer :: i

    ebind_per_nucleon(jn)   = 0.782347d0
    ebind_per_nucleon(jp)   = 0.0d0
    ebind_per_nucleon(jhe4)   = 7.073915d0
    ebind_per_nucleon(jc12)   = 7.680144d0
    ebind_per_nucleon(jne20)   = 8.032240d0
    ebind_per_nucleon(jna23)   = 8.111493d0
    ebind_per_nucleon(jmg23)   = 7.901104d0
    
    aion(jn)   = 1.000000d+00
    aion(jp)   = 1.000000d+00
    aion(jhe4)   = 4.000000d+00
    aion(jc12)   = 1.200000d+01
    aion(jne20)   = 2.000000d+01
    aion(jna23)   = 2.300000d+01
    aion(jmg23)   = 2.300000d+01

    zion(jn)   =     0.000000d+00
    zion(jp)   =     1.000000d+00
    zion(jhe4)   =   2.000000d+00
    zion(jc12)   =   6.000000d+00
    zion(jne20)   = 10.00000d+00
    zion(jna23)   = 11.00000d+00
    zion(jmg23)   = 12.00000d+00
    
    do i = 1, nspec
      ebind(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do
  end subroutine init_network

end module network
