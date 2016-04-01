module net_rates
  use physical_constants
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network, only: nspec, aion, zion, net_meta
  use table_rates

  implicit none

  integer, parameter :: nreact = 4
  logical :: screen_reaclib = .true.

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_1
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_2
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_3
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_4

  type :: reaclib_coefs
    ! how many rates in this reaction
    integer :: m
    ! coefficient array
    ! ! First Index: Rate coefficients 1..7
    ! ! Second Index: Reaction index 1..m
    double precision, dimension(:,:), allocatable :: ctemp
  end type reaclib_coefs

  ! Reaction multiplicities (how many rates contribute)
  ! array indexed by the Reactions indices above
  integer, dimension(nreact) :: rate_mult = (/ &
       1, &
       1, &
       1, &
       1 /)

  ! Should these reactions be screened?
  logical, dimension(nreact) :: do_screening = (/ &
       .true., &
       .true., &
       .true., &
       .false. /)

contains

  subroutine init_reaclib_pars(pfile_unit)
    integer, intent(in) :: pfile_unit
    
    namelist /reaclibpars/ screen_reaclib
    
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=reaclibpars)
  end subroutine init_reaclib_pars
  
  subroutine init_reaclib()
    
    allocate( ctemp_rate_1(7, rate_mult(1)) )
    ! c12_c12a_ne20
    ctemp_rate_1(:, 1) = (/  &
        6.128630d+01, &
        0.000000d+00, &
        -8.416500d+01, &
        -1.566270d+00, &
        -7.360840d-02, &
        -7.279700d-02, &
        -6.666670d-01 /)

    allocate( ctemp_rate_2(7, rate_mult(2)) )
    ! c12_c12n_mg23
    ctemp_rate_2(:, 1) = (/  &
        -1.280560d+01, &
        -3.014850d+01, &
        0.000000d+00, &
        1.148260d+01, &
        1.828490d+00, &
        -3.484400d-01, &
        0.000000d+00 /)

    allocate( ctemp_rate_3(7, rate_mult(3)) )
    ! c12_c12p_na23
    ctemp_rate_3(:, 1) = (/  &
        6.096490d+01, &
        0.000000d+00, &
        -8.416500d+01, &
        -1.419100d+00, &
        -1.146190d-01, &
        -7.030700d-02, &
        -6.666670d-01 /)

    allocate( ctemp_rate_4(7, rate_mult(4)) )
    ! n_p
    ctemp_rate_4(:, 1) = (/  &
        -6.781610d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00, &
        0.000000d+00 /)

  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate_1 )
    deallocate( ctemp_rate_2 )
    deallocate( ctemp_rate_3 )
    deallocate( ctemp_rate_4 )
  end subroutine term_reaclib

  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init
    !    <add_screening>(2)

    call add_screening_factor(zion(net_meta%ic12), aion(net_meta%ic12), &
         zion(net_meta%ic12), aion(net_meta%ic12))

    call add_screening_factor(zion(net_meta%ic12), aion(net_meta%ic12), &
         zion(net_meta%ic12), aion(net_meta%ic12))

    call add_screening_factor(zion(net_meta%ic12), aion(net_meta%ic12), &
         zion(net_meta%ic12), aion(net_meta%ic12))

!    call add_screening_factor(zion(net_meta%in), aion(net_meta%in), &
!         zion(net_meta%ip), aion(net_meta%ip)) ! Is this right for n --> p ?
        
    call screening_init()    
  end subroutine net_screening_init

  subroutine rate_evaluate(state, rhoy, temp, iwhich, ratevec)
    type(plasma_state), intent(in) :: state
    double precision, intent(in) :: temp, rhoy
    integer, intent(in) :: iwhich

    double precision, intent(out), dimension(:) :: ratevec
    ! ratevec(1) = rate     , the reaction rate
    ! ratevec(2) = drate_dt , the Temperature derivative of rate
    ! ratevec(3) = scor     , the screening factor
    ! ratevec(4) = dscor_dt , the Temperature derivative of scor

    double precision  :: rate, scor ! Rate and Screening Factor
    double precision  :: drate_dt, dscor_dt ! Temperature derivatives
    double precision :: dscor_dd
    double precision, pointer :: ctemp(:,:)
    double precision :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m
    
    ri = 0.0d0
    rate = 0.0d0
    drate_dt = 0.0d0
    irate = 0.0d0
    dirate_dt = 0.0d0
    T9 = temp/1.0d9
    T9_exp = 0.0d0
    scor = 1.0d0
    dscor_dt = 0.0d0
    dscor_dd = 0.0d0

    if (iwhich == 1) then
      ctemp => ctemp_rate_1
    else if (iwhich == 2) then
      ctemp => ctemp_rate_2
    else if (iwhich == 3) then
      ctemp => ctemp_rate_3
    else if (iwhich == 4) then
      ctemp => ctemp_rate_4
    end if

    m = rate_mult(iwhich)
    do i = 1, m
      lnirate = ctemp(1,i) + ctemp(7,i) * LOG(T9)
      dlnirate_dt = ctemp(7,i)/T9
      do j = 2, 6
        T9_exp = (2.0d0*dble(j-1)-5.0d0)/3.0d0 
        lnirate = lnirate + ctemp(j,i) * T9**T9_exp
        dlnirate_dt = dlnirate_dt + T9_exp * ctemp(j,i) * T9**(T9_exp-1.0d0)
      end do
      irate = EXP(lnirate)
      rate = rate + irate
      dirate_dt = irate * dlnirate_dt/1.0d9
      drate_dt = drate_dt + dirate_dt
    end do
    
    if ( screen_reaclib .and. do_screening(iwhich) ) then
      call screen5(state, iwhich, scor, dscor_dt, dscor_dd)
    end if

    ratevec(1) = rate
    ratevec(2) = drate_dt
    ratevec(3) = scor
    ratevec(4) = dscor_dt
  end subroutine rate_evaluate

end module net_rates
