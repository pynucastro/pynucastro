module net_rates
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network

  implicit none

  logical, parameter :: screen_reaclib = .true.
  
  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start  
  integer, dimension(nrat_reaclib) :: rate_start_idx = (/ &
    1, &
    3, &
    4, &
    7, &
    11, &
    12, &
    16, &
    18 /)
  
  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, dimension(nrat_reaclib) :: rate_extra_mult = (/ &
    1, &
    0, &
    2, &
    3, &
    0, &
    3, &
    1, &
    0 /)

  ! Should these reactions be screened?
  logical, dimension(nrat_reaclib) :: do_screening = (/ &
    .true., &
    .false., &
    .true., &
    .true., &
    .false., &
    .true., &
    .true., &
    .false. /)

contains

  subroutine init_reaclib()

    allocate( ctemp_rate(7, 18) )
    ! c12_pg_n13
    ctemp_rate(:, 1) = (/  &
        1.71482000000000d+01, &
        0.00000000000000d+00, &
        -1.36920000000000d+01, &
        -2.30881000000000d-01, &
        4.44362000000000d+00, &
        -3.15898000000000d+00, &
        -6.66667000000000d-01 /)

    ctemp_rate(:, 2) = (/  &
        1.75428000000000d+01, &
        -3.77849000000000d+00, &
        -5.10735000000000d+00, &
        -2.24111000000000d+00, &
        1.48883000000000d-01, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 /)

    ! n13_c13
    ctemp_rate(:, 3) = (/  &
        -6.76010000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 /)

    ! c13_pg_n14
    ctemp_rate(:, 4) = (/  &
        1.85155000000000d+01, &
        0.00000000000000d+00, &
        -1.37200000000000d+01, &
        -4.50018000000000d-01, &
        3.70823000000000d+00, &
        -1.70545000000000d+00, &
        -6.66667000000000d-01 /)

    ctemp_rate(:, 5) = (/  &
        1.39637000000000d+01, &
        -5.78147000000000d+00, &
        0.00000000000000d+00, &
        -1.96703000000000d-01, &
        1.42126000000000d-01, &
        -2.38912000000000d-02, &
        -1.50000000000000d+00 /)

    ctemp_rate(:, 6) = (/  &
        1.51825000000000d+01, &
        -1.35543000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 /)

    ! n14_pg_o15
    ctemp_rate(:, 7) = (/  &
        1.70100000000000d+01, &
        0.00000000000000d+00, &
        -1.51930000000000d+01, &
        -1.61954000000000d-01, &
        -7.52123000000000d+00, &
        -9.87565000000000d-01, &
        -6.66667000000000d-01 /)

    ctemp_rate(:, 8) = (/  &
        6.73578000000000d+00, &
        -4.89100000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        6.82000000000000d-02 /)

    ctemp_rate(:, 9) = (/  &
        7.65444000000000d+00, &
        -2.99800000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 /)

    ctemp_rate(:, 10) = (/  &
        2.01169000000000d+01, &
        0.00000000000000d+00, &
        -1.51930000000000d+01, &
        -4.63975000000000d+00, &
        9.73458000000000d+00, &
        -9.55051000000000d+00, &
        3.33333000000000d-01 /)

    ! o15_n15
    ctemp_rate(:, 11) = (/  &
        -5.17053000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 /)

    ! n15_pa_c12
    ctemp_rate(:, 12) = (/  &
        2.74764000000000d+01, &
        0.00000000000000d+00, &
        -1.52530000000000d+01, &
        1.59318000000000d+00, &
        2.44790000000000d+00, &
        -2.19708000000000d+00, &
        -6.66667000000000d-01 /)

    ctemp_rate(:, 13) = (/  &
        -6.57522000000000d+00, &
        -1.16380000000000d+00, &
        0.00000000000000d+00, &
        2.27105000000000d+01, &
        -2.90707000000000d+00, &
        2.05754000000000d-01, &
        -1.50000000000000d+00 /)

    ctemp_rate(:, 14) = (/  &
        2.08972000000000d+01, &
        -7.40600000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 /)

    ctemp_rate(:, 15) = (/  &
        -4.87347000000000d+00, &
        -2.02117000000000d+00, &
        0.00000000000000d+00, &
        3.08497000000000d+01, &
        -8.50433000000000d+00, &
        -1.54426000000000d+00, &
        -1.50000000000000d+00 /)

    ! n13_pg_o14
    ctemp_rate(:, 16) = (/  &
        1.81356000000000d+01, &
        0.00000000000000d+00, &
        -1.51676000000000d+01, &
        9.55166000000000d-02, &
        3.06590000000000d+00, &
        -5.07339000000000d-01, &
        -6.66667000000000d-01 /)

    ctemp_rate(:, 17) = (/  &
        1.09971000000000d+01, &
        -6.12602000000000d+00, &
        1.57122000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 /)

    ! o14_n14
    ctemp_rate(:, 18) = (/  &
        -4.62354000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 /)



  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate )
  end subroutine term_reaclib

  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn15), aion(jn15))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn13), aion(jn13))


    call screening_init()    
  end subroutine net_screening_init

  subroutine reaclib_evaluate(pstate, temp, iwhich, reactvec)

    implicit none
    
    type(plasma_state), intent(in) :: pstate
    double precision, intent(in) :: temp
    integer, intent(in) :: iwhich

    double precision, intent(inout) :: reactvec(6)
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
    double precision :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m, istart

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
    m = rate_extra_mult(iwhich)

    istart = rate_start_idx(iwhich)

    do i = 0, m
       lnirate = ctemp_rate(1, istart+i) + ctemp_rate(7, istart+i) * LOG(T9)
       dlnirate_dt = ctemp_rate(7, istart+i)/T9
       do j = 2, 6
          T9_exp = (2.0d0*dble(j-1)-5.0d0)/3.0d0 
          lnirate = lnirate + ctemp_rate(j, istart+i) * T9**T9_exp
          dlnirate_dt = dlnirate_dt + &
               T9_exp * ctemp_rate(j, istart+i) * T9**(T9_exp-1.0d0)
       end do
       ! If the rate will be in the approx. interval [0.0, 1.0E-100], replace by 0.0
       ! This avoids issues with passing very large negative values to EXP
       ! and getting results between 0.0 and 1.0E-308, the limit for IEEE 754.
       ! And avoids SIGFPE in CVODE due to tiny rates.
       lnirate = max(lnirate, -230.0d0)
       irate = EXP(lnirate)
       rate = rate + irate
       dirate_dt = irate * dlnirate_dt/1.0d9
       drate_dt = drate_dt + dirate_dt
    end do

    if ( screen_reaclib .and. do_screening(iwhich) ) then
       call screen5(pstate, iwhich, scor, dscor_dt, dscor_dd)
    end if

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = scor
    reactvec(i_dscor_dt) = dscor_dt
    reactvec(i_dqweak)   = 0.0d0
    reactvec(i_epart)    = 0.0d0

    ! write(*,*) '----------------------------------------'
    ! write(*,*) 'IWHICH: ', iwhich
    ! write(*,*) 'reactvec(i_rate)', reactvec(i_rate)
    ! write(*,*) 'reactvec(i_drate_dt)', reactvec(i_drate_dt)
    ! write(*,*) 'reactvec(i_scor)', reactvec(i_scor)    
    ! write(*,*) 'reactvec(i_dscor_dt)', reactvec(i_dscor_dt)
    ! write(*,*) 'reactvec(i_dqweak)', reactvec(i_dqweak)
    ! write(*,*) 'reactvec(i_epart)', reactvec(i_epart)
    ! write(*,*) '----------------------------------------'

  end subroutine reaclib_evaluate
  
end module net_rates
