module net_rates
  use physical_constants
  use table_rates

  implicit none

  integer, parameter :: number_equations = 9
  integer, parameter :: number_nuclides = 8
  integer, parameter :: number_reactions = 6

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_1
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_2
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_3
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_6

  type :: net_initials_t
    ! Initial abundances
    double precision :: yn
    double precision :: yp
    double precision :: yhe4
    double precision :: yc12
    double precision :: yne20
    double precision :: yne23
    double precision :: yna23
    double precision :: ymg23
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
    integer :: in   = 1
    integer :: ip   = 2
    integer :: ihe4   = 3
    integer :: ic12   = 4
    integer :: ine20   = 5
    integer :: ine23   = 6
    integer :: ina23   = 7
    integer :: img23   = 8
    ! Energy Generation Rate
    integer :: ienuc   = 9

    ! Reactions
    integer :: k_c12_c12a_ne20   = 1
    integer :: k_c12_c12n_mg23   = 2
    integer :: k_c12_c12p_na23   = 3
    integer :: k_na23_ne23   = 4
    integer :: k_ne23_na23   = 5
    integer :: k_n_p   = 6

    ! Reaction multiplicities (how many rates contribute)
    ! array indexed by the Reactions indices above
    integer, dimension(number_reactions) :: rate_mult = (/ &
    1, &
    1, &
    1, &
    0, &
    0, &
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

    self%ebind_per_nucleon(self%in)   = 0.0d0
    self%ebind_per_nucleon(self%ip)   = 0.0d0
    self%ebind_per_nucleon(self%ihe4)   = 0.0d0
    self%ebind_per_nucleon(self%ic12)   = 0.0d0
    self%ebind_per_nucleon(self%ine20)   = 0.0d0
    self%ebind_per_nucleon(self%ine23)   = 0.0d0
    self%ebind_per_nucleon(self%ina23)   = 0.0d0
    self%ebind_per_nucleon(self%img23)   = 0.0d0

    self%anuc(self%in)   = 1.000000d+00
    self%anuc(self%ip)   = 1.000000d+00
    self%anuc(self%ihe4)   = 4.000000d+00
    self%anuc(self%ic12)   = 1.200000d+01
    self%anuc(self%ine20)   = 2.000000d+01
    self%anuc(self%ine23)   = 2.300000d+01
    self%anuc(self%ina23)   = 2.300000d+01
    self%anuc(self%img23)   = 2.300000d+01

    do i = 1, self%nnuc
      self%ebind(i) = self%ebind_per_nucleon(i) * self%anuc(i) * pc%erg_per_MeV
    end do

    allocate( ctemp_rate_1(7, self%rate_mult(1)) )
    ! c12_c12a_ne20
    ctemp_rate_1(:, 1) = (/  &
        6.128630d+01, &
        0.000000d+00, &
        -8.416500d+01, &
        -1.566270d+00, &
        -7.360840d-02, &
        -7.279700d-02, &
        -6.666670d-01 /)

    allocate( ctemp_rate_2(7, self%rate_mult(2)) )
    ! c12_c12n_mg23
    ctemp_rate_2(:, 1) = (/  &
        -1.280560d+01, &
        -3.014850d+01, &
        0.000000d+00, &
        1.148260d+01, &
        1.828490d+00, &
        -3.484400d-01, &
        0.000000d+00 /)

    allocate( ctemp_rate_3(7, self%rate_mult(3)) )
    ! c12_c12p_na23
    ctemp_rate_3(:, 1) = (/  &
        6.096490d+01, &
        0.000000d+00, &
        -8.416500d+01, &
        -1.419100d+00, &
        -1.146190d-01, &
        -7.030700d-02, &
        -6.666670d-01 /)

    allocate( ctemp_rate_6(7, self%rate_mult(6)) )
    ! n_p
    ctemp_rate_6(:, 1) = (/  &
        -6.781610d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00 /)

    call init_table_meta()




  end subroutine init_net_info

  subroutine term_net_info(self)
    class(net_info) :: self
    deallocate( ctemp_rate_1 )
    deallocate( ctemp_rate_2 )
    deallocate( ctemp_rate_3 )
    deallocate( ctemp_rate_6 )
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
      call table_meta(j_na23_ne23)%bl_lookup(rhoy, temp, jtab_rate, rate)
      return
    else if (iwhich == 5) then
      call table_meta(j_ne23_na23)%bl_lookup(rhoy, temp, jtab_rate, rate)
      return
    else if (iwhich == 6) then
      ctemp => ctemp_rate_6
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