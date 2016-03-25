module net_rates
  use physical_constants
  use table_rates

  implicit none

  integer, parameter :: number_equations = 10
  integer, parameter :: number_nuclides = 9
  integer, parameter :: number_reactions = 8

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_1
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_2
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_3
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_4
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_5
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_6
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_7
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_8

  type :: net_initials_t
    ! Initial abundances
    double precision :: yp
    double precision :: yhe4
    double precision :: yc12
    double precision :: yc13
    double precision :: yn13
    double precision :: yn14
    double precision :: yn15
    double precision :: yo14
    double precision :: yo15
    double precision :: yenuc = 0.0d0
  end type net_initials_t

  type :: reaclib_coefs
    ! how many rates in this reaction
    integer :: m
    ! coefficient array
    ! ! First Index: Rate coefficients 1..7
    ! ! Second Index: Reaction index 1..m
    double precision, dimension(:,:), allocatable :: ctemp
  end type reaclib_coefs

  type :: net_info
    integer :: neqs  = number_equations
    integer :: nnuc  = number_nuclides
    integer :: nrxn  = number_reactions
    ! Nuclides
    integer :: ip   = 1
    integer :: ihe4   = 2
    integer :: ic12   = 3
    integer :: ic13   = 4
    integer :: in13   = 5
    integer :: in14   = 6
    integer :: in15   = 7
    integer :: io14   = 8
    integer :: io15   = 9
    ! Energy Generation Rate
    integer :: ienuc   = 10

    ! Reactions
    integer :: k_c12_pg_n13   = 1
    integer :: k_n13_c13   = 2
    integer :: k_c13_pg_n14   = 3
    integer :: k_n14_pg_o15   = 4
    integer :: k_o15_n15   = 5
    integer :: k_n15_pa_c12   = 6
    integer :: k_n13_pg_o14   = 7
    integer :: k_o14_n14   = 8

    ! Reaction multiplicities (how many rates contribute)
    ! array indexed by the Reactions indices above
    integer, dimension(number_reactions) :: rate_mult = (/ &
    2, &
    1, &
    3, &
    4, &
    1, &
    4, &
    2, &
    1 /)

    ! Binding Energies Per Nucleon (MeV)
    double precision, dimension(number_nuclides) :: ebind_per_nucleon

    ! Nucleon number A
    double precision, dimension(number_nuclides) :: anuc

    ! Binding Energies (ergs)
    double precision, dimension(number_nuclides) :: ebind

  contains
    procedure :: initialize => init_net_info
    procedure :: terminate => term_net_info
    procedure :: evaluate => rates_eval
  end type net_info

  type(net_info), target, save :: net_meta
  type(net_initials_t), save :: net_initial_abundances

contains

  subroutine init_net_pars(pfile_unit)
    integer, intent(in) :: pfile_unit

    namelist /netpars/ net_initial_abundances
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=netpars)
  end subroutine init_net_pars

  subroutine init_net_info(self)
    class(net_info) :: self
    type(reaclib_coefs), pointer :: prate
    integer :: i

    self%ebind_per_nucleon(self%ip)   = 0.0d0
    self%ebind_per_nucleon(self%ihe4)   = 0.0d0
    self%ebind_per_nucleon(self%ic12)   = 0.0d0
    self%ebind_per_nucleon(self%ic13)   = 0.0d0
    self%ebind_per_nucleon(self%in13)   = 0.0d0
    self%ebind_per_nucleon(self%in14)   = 0.0d0
    self%ebind_per_nucleon(self%in15)   = 0.0d0
    self%ebind_per_nucleon(self%io14)   = 0.0d0
    self%ebind_per_nucleon(self%io15)   = 0.0d0

    self%anuc(self%ip)   = 1.000000d+00
    self%anuc(self%ihe4)   = 4.000000d+00
    self%anuc(self%ic12)   = 1.200000d+01
    self%anuc(self%ic13)   = 1.300000d+01
    self%anuc(self%in13)   = 1.300000d+01
    self%anuc(self%in14)   = 1.400000d+01
    self%anuc(self%in15)   = 1.500000d+01
    self%anuc(self%io14)   = 1.400000d+01
    self%anuc(self%io15)   = 1.500000d+01

    do i = 1, self%nnuc
      self%ebind(i) = self%ebind_per_nucleon(i) * self%anuc(i) * pc%erg_per_MeV
    end do

    allocate( ctemp_rate_1(7, self%rate_mult(1)) )
    ! c12_pg_n13
    ctemp_rate_1(:, 1) = (/  &
        1.714820d+01, &
        0.000000d+00, &
        -1.369200d+01, &
        -2.308810d-01, &
        4.443620d+00, &
        -3.158980d+00, &
        -6.666670d-01 /)

    ctemp_rate_1(:, 2) = (/  &
        1.754280d+01, &
        -3.778490d+00, &
        -5.107350d+00, &
        -2.241110d+00, &
        1.488830d-01, &
        0.000000d+00, &
        -1.500000d+00 /)

    allocate( ctemp_rate_2(7, self%rate_mult(2)) )
    ! n13_c13
    ctemp_rate_2(:, 1) = (/  &
        -6.760100d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00 /)

    allocate( ctemp_rate_3(7, self%rate_mult(3)) )
    ! c13_pg_n14
    ctemp_rate_3(:, 1) = (/  &
        1.851550d+01, &
        0.000000d+00, &
        -1.372000d+01, &
        -4.500180d-01, &
        3.708230d+00, &
        -1.705450d+00, &
        -6.666670d-01 /)

    ctemp_rate_3(:, 2) = (/  &
        1.396370d+01, &
        -5.781470d+00, &
        0.000000d+00, &
        -1.967030d-01, &
        1.421260d-01, &
        -2.389120d-02, &
        -1.500000d+00 /)

    ctemp_rate_3(:, 3) = (/  &
        1.518250d+01, &
        -1.355430d+01, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        -1.500000d+00 /)

    allocate( ctemp_rate_4(7, self%rate_mult(4)) )
    ! n14_pg_o15
    ctemp_rate_4(:, 1) = (/  &
        1.701000d+01, &
        0.000000d+00, &
        -1.519300d+01, &
        -1.619540d-01, &
        -7.521230d+00, &
        -9.875650d-01, &
        -6.666670d-01 /)

    ctemp_rate_4(:, 2) = (/  &
        6.735780d+00, &
        -4.891000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        6.820000d-02 /)

    ctemp_rate_4(:, 3) = (/  &
        7.654440d+00, &
        -2.998000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        -1.500000d+00 /)

    ctemp_rate_4(:, 4) = (/  &
        2.011690d+01, &
        0.000000d+00, &
        -1.519300d+01, &
        -4.639750d+00, &
        9.734580d+00, &
        -9.550510d+00, &
        3.333330d-01 /)

    allocate( ctemp_rate_5(7, self%rate_mult(5)) )
    ! o15_n15
    ctemp_rate_5(:, 1) = (/  &
        -5.170530d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00 /)

    allocate( ctemp_rate_6(7, self%rate_mult(6)) )
    ! n15_pa_c12
    ctemp_rate_6(:, 1) = (/  &
        2.747640d+01, &
        0.000000d+00, &
        -1.525300d+01, &
        1.593180d+00, &
        2.447900d+00, &
        -2.197080d+00, &
        -6.666670d-01 /)

    ctemp_rate_6(:, 2) = (/  &
        -6.575220d+00, &
        -1.163800d+00, &
        0.000000d+00, &
        2.271050d+01, &
        -2.907070d+00, &
        2.057540d-01, &
        -1.500000d+00 /)

    ctemp_rate_6(:, 3) = (/  &
        2.089720d+01, &
        -7.406000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        -1.500000d+00 /)

    ctemp_rate_6(:, 4) = (/  &
        -4.873470d+00, &
        -2.021170d+00, &
        0.000000d+00, &
        3.084970d+01, &
        -8.504330d+00, &
        -1.544260d+00, &
        -1.500000d+00 /)

    allocate( ctemp_rate_7(7, self%rate_mult(7)) )
    ! n13_pg_o14
    ctemp_rate_7(:, 1) = (/  &
        1.813560d+01, &
        0.000000d+00, &
        -1.516760d+01, &
        9.551660d-02, &
        3.065900d+00, &
        -5.073390d-01, &
        -6.666670d-01 /)

    ctemp_rate_7(:, 2) = (/  &
        1.099710d+01, &
        -6.126020d+00, &
        1.571220d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        -1.500000d+00 /)

    allocate( ctemp_rate_8(7, self%rate_mult(8)) )
    ! o14_n14
    ctemp_rate_8(:, 1) = (/  &
        -4.623540d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00 /)





  end subroutine init_net_info

  subroutine term_net_info(self)
    class(net_info) :: self
    deallocate( ctemp_rate_1 )
    deallocate( ctemp_rate_2 )
    deallocate( ctemp_rate_3 )
    deallocate( ctemp_rate_4 )
    deallocate( ctemp_rate_5 )
    deallocate( ctemp_rate_6 )
    deallocate( ctemp_rate_7 )
    deallocate( ctemp_rate_8 )
  end subroutine term_net_info

  subroutine rates_eval(self, rhoy, temp, iwhich, rate)
    class(net_info) :: self
    double precision, intent(in) :: temp, rhoy
    integer, intent(in) :: iwhich
    double precision, intent(out) :: rate
    double precision, pointer :: ctemp(:,:)

    double precision :: ri, T9, lnrate
    integer :: i, j, m

    ri = 0.0d0
    rate = 0.0d0
    T9 = temp/1.0d9

    if (iwhich == 1) then
      ctemp => ctemp_rate_1
    else if (iwhich == 2) then
      ctemp => ctemp_rate_2
    else if (iwhich == 3) then
      ctemp => ctemp_rate_3
    else if (iwhich == 4) then
      ctemp => ctemp_rate_4
    else if (iwhich == 5) then
      ctemp => ctemp_rate_5
    else if (iwhich == 6) then
      ctemp => ctemp_rate_6
    else if (iwhich == 7) then
      ctemp => ctemp_rate_7
    else if (iwhich == 8) then
      ctemp => ctemp_rate_8
    end if

    m = self%rate_mult(iwhich)
    do i = 1, m
      lnrate = ctemp(1,i) + ctemp(7,i) * LOG(T9)
      do j = 2, 6
        lnrate = lnrate + ctemp(j,i) * T9**((2.0d0*dble(j-1)-5.0d0)/3.0d0)
      end do
      rate = rate + EXP(lnrate)
    end do
  end subroutine rates_eval

end module net_rates