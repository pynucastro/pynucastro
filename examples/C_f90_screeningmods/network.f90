module network
  use physical_constants, only: pc
  
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
  
  type :: net_initials_t
    ! Initial abundances (Y = X/A)
    double precision :: yn
    double precision :: yp
    double precision :: yhe4
    double precision :: yc12
    double precision :: yne20
    double precision :: yna23
    double precision :: ymg23
    double precision :: yenuc = 0.0d0
  end type net_initials_t

  type :: net_info
    ! Nuclides
    integer :: in   = 1
    integer :: ip   = 2
    integer :: ihe4   = 3
    integer :: ic12   = 4
    integer :: ine20   = 5
    integer :: ina23   = 6
    integer :: img23   = 7
    ! Energy Generation Rate
    integer :: ienuc   = 8

    ! Reactions
    integer :: k_c12_c12a_ne20   = 1
    integer :: k_c12_c12n_mg23   = 2
    integer :: k_c12_c12p_na23   = 3
    integer :: k_n_p   = 4
  end type net_info

  type(net_info) :: net_meta
  type(net_initials_t), save :: net_initial_abundances
  
contains
  
  subroutine init_net_pars(pfile_unit)
    integer, intent(in) :: pfile_unit

    namelist /netpars/ net_initial_abundances

    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=netpars)
  end subroutine init_net_pars

  subroutine init_network()
    integer :: i

    ebind_per_nucleon(net_meta%in)   = 0.782347d0
    ebind_per_nucleon(net_meta%ip)   = 0.0d0
    ebind_per_nucleon(net_meta%ihe4)   = 7.073915d0
    ebind_per_nucleon(net_meta%ic12)   = 7.680144d0
    ebind_per_nucleon(net_meta%ine20)   = 8.032240d0
    ebind_per_nucleon(net_meta%ina23)   = 8.111493d0
    ebind_per_nucleon(net_meta%img23)   = 7.901104d0
    
    aion(net_meta%in)   = 1.000000d+00
    aion(net_meta%ip)   = 1.000000d+00
    aion(net_meta%ihe4)   = 4.000000d+00
    aion(net_meta%ic12)   = 1.200000d+01
    aion(net_meta%ine20)   = 2.000000d+01
    aion(net_meta%ina23)   = 2.300000d+01
    aion(net_meta%img23)   = 2.300000d+01

    zion(net_meta%in)   =     0.000000d+00
    zion(net_meta%ip)   =     1.000000d+00
    zion(net_meta%ihe4)   =   2.000000d+00
    zion(net_meta%ic12)   =   6.000000d+00
    zion(net_meta%ine20)   = 10.00000d+00
    zion(net_meta%ina23)   = 11.00000d+00
    zion(net_meta%img23)   = 12.00000d+00
    
    do i = 1, nspec
      ebind(i) = ebind_per_nucleon(i) * aion(i) * pc%erg_per_MeV
    end do
  end subroutine init_network

end module network
