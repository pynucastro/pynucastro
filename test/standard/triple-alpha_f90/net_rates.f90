module net_rates
  use physical_constants
  use table_rates

  implicit none

  integer, parameter :: number_equations = 3
  integer, parameter :: number_nuclides = 2
  integer, parameter :: number_reactions = 2

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_1
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_2

  type :: net_initials_t
    ! Initial abundances
    double precision :: yhe4
    double precision :: yc12
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
    integer :: ihe4   = 1
    integer :: ic12   = 2
    ! Energy Generation Rate
    integer :: ienuc   = 3

    ! Reactions
    integer :: k_c12_gaa_he4   = 1
    integer :: k_he4_aag_c12   = 2

    ! Reaction multiplicities (how many rates contribute)
    ! array indexed by the Reactions indices above
    integer, dimension(number_reactions) :: rate_mult = (/ &
    3, &
    3 /)

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

    self%ebind_per_nucleon(self%ihe4)   = 0.0d0
    self%ebind_per_nucleon(self%ic12)   = 0.0d0

    self%anuc(self%ihe4)   = 4.000000d+00
    self%anuc(self%ic12)   = 1.200000d+01

    do i = 1, self%nnuc
      self%ebind(i) = self%ebind_per_nucleon(i) * self%anuc(i) * pc%erg_per_MeV
    end do

    allocate( ctemp_rate_1(7, self%rate_mult(1)) )
    ! c12_gaa_he4
    ctemp_rate_1(:, 1) = (/  &
        3.495610d+01, &
        -8.544720d+01, &
        -2.357000d+01, &
        2.048860d+01, &
        -1.298820d+01, &
        -2.000000d+01, &
        8.333300d-01 /)

    ctemp_rate_1(:, 2) = (/  &
        4.577340d+01, &
        -8.442270d+01, &
        -3.706000d+01, &
        2.934930d+01, &
        -1.155070d+02, &
        -1.000000d+01, &
        1.666670d+00 /)

    ctemp_rate_1(:, 3) = (/  &
        2.239400d+01, &
        -8.854930d+01, &
        -1.349000d+01, &
        2.142590d+01, &
        -1.347690d+00, &
        8.798160d-02, &
        -1.016530d+01 /)

    allocate( ctemp_rate_2(7, self%rate_mult(2)) )
    ! he4_aag_c12
    ctemp_rate_2(:, 1) = (/  &
        -9.710520d-01, &
        0.000000d+00, &
        -3.706000d+01, &
        2.934930d+01, &
        -1.155070d+02, &
        -1.000000d+01, &
        -1.333330d+00 /)

    ctemp_rate_2(:, 2) = (/  &
        -2.435050d+01, &
        -4.126560d+00, &
        -1.349000d+01, &
        2.142590d+01, &
        -1.347690d+00, &
        8.798160d-02, &
        -1.316530d+01 /)

    ctemp_rate_2(:, 3) = (/  &
        -1.178840d+01, &
        -1.024460d+00, &
        -2.357000d+01, &
        2.048860d+01, &
        -1.298820d+01, &
        -2.000000d+01, &
        -2.166670d+00 /)





  end subroutine init_net_info

  subroutine term_net_info(self)
    class(net_info) :: self
    deallocate( ctemp_rate_1 )
    deallocate( ctemp_rate_2 )
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