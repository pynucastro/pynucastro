module net_rates
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network
  use table_rates

  implicit none

  integer, parameter :: nreact = 6
  logical :: screen_reaclib = .true.
  
  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_1
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_2
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_3
  double precision, target, dimension(:,:), allocatable :: ctemp_rate_4

  ! type :: reaclib_coefs
  !   ! how many rates in this reaction
  !   integer :: m
  !   ! coefficient array
  !   ! ! First Index: Rate coefficients 1..7
  !   ! ! Second Index: Reaction index 1..m
  !   double precision, dimension(:,:), allocatable :: ctemp
  ! end type reaclib_coefs

  type :: ctemp_ptr
     double precision, dimension(:,:), pointer :: p
  end type ctemp_ptr
  
  ! Declare an array of pointers to ctemp arrays
  type(ctemp_ptr), dimension(4) :: ctemp_point
  
  ! Reaction multiplicities (how many rates contribute)
  ! Array indexed by the Reactions indices above.
  ! If the entry is negative, interpret as an index into table_meta
  ! and calculate the rate using module table_rates.
  integer, dimension(nreact) :: rate_mult = (/ &
    1, &
    1, &
    1, &
    1, &
    -1, &
    -2 /)

  ! Should these reactions be screened?
  logical, dimension(nreact) :: do_screening = (/ &
    .true., &
    .true., &
    .true., &
    .false., &
    .false., &
    .false. /)

contains

  subroutine init_reaclib_pars(pfile_unit)
    integer, intent(in) :: pfile_unit
    
    namelist /reaclibpars/ screen_reaclib
    
    rewind(unit=pfile_unit)
    read(unit=pfile_unit, nml=reaclibpars)
  end subroutine init_reaclib_pars

  subroutine init_reaclib()
    ctemp_point(1)%p => ctemp_rate_1    
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

    ctemp_point(2)%p => ctemp_rate_2
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

    ctemp_point(3)%p => ctemp_rate_3
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
    
    ctemp_point(4)%p => ctemp_rate_4
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

    call init_table_meta()

  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate_1 )
    deallocate( ctemp_rate_2 )
    deallocate( ctemp_rate_3 )
    deallocate( ctemp_rate_4 )
  end subroutine term_reaclib

  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call screening_init()    
  end subroutine net_screening_init

  subroutine rate_evaluate(state, rhoy, temp, iwhich, reactvec)
    type(plasma_state), intent(in) :: state
    double precision, intent(in) :: temp, rhoy
    integer, intent(in) :: iwhich

    double precision, intent(out), dimension(:) :: reactvec
    ! reactvec(1) = rate     , the reaction rate
    ! reactvec(2) = drate_dt , the Temperature derivative of rate
    ! reactvec(3) = scor     , the screening factor
    ! reactvec(4) = dscor_dt , the Temperature derivative of scor
    ! reactvec(5) = dqweak   , the weak reaction dq-value (ergs)
    !                          (This accounts for modification of the reaction Q
    !                           due to the local density and temperature of the plasma.
    !                           For Reaclib rates, this is 0.0d0.)
    ! reactvec(6) = epart    , the particle energy generation rate (ergs/s)
    ! NOTE: The particle energy generation rate (returned in ergs/s)
    !       is the contribution to enuc from non-ion particles associated
    !       with the reaction.
    !       For example, this accounts for neutrino energy losses
    !       in weak reactions and/or gamma heating of the plasma
    !       from nuclear transitions in daughter nuclei.

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

    ! Use reaction multiplicities to tell whether the rate is Reaclib
    m = rate_mult(iwhich)
    if ( m .gt. 0 ) then
       ! This must be a Reaclib rate
       ctemp = ctemp_point(iwhich)%p
       
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
 
       reactvec(1) = rate
       reactvec(2) = drate_dt
       reactvec(3) = scor
       reactvec(4) = dscor_dt
       reactvec(5) = 0.0d0
       reactvec(6) = 0.0d0
    else
       ! This is a Table rate (unless m=0)
       ! table_rates returns dq in reactvec(5)
       ! reflecting the dependence of Q on dens*ye and temperature.
       call table_meta(-m)%reaction(rhoy, temp, reactvec)
    end if
  end subroutine rate_evaluate

end module net_rates
