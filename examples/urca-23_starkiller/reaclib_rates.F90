module reaclib_rates

  use amrex_fort_module, only: rt => amrex_real
  use screening_module, only: add_screening_factor, &
                              screening_init, screening_finalize, &
                              plasma_state, fill_plasma_state
  use network

  implicit none

  logical, parameter :: screen_reaclib = .true.

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  real(rt), allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start
  integer, allocatable :: rate_start_idx(:)

  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, allocatable :: rate_extra_mult(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ctemp_rate, rate_start_idx, rate_extra_mult
#endif

  !$acc declare create(ctemp_rate, rate_start_idx, rate_extra_mult)
  !$acc declare copyin(screen_reaclib)

contains

  subroutine init_reaclib()

    implicit none

    integer :: unit, ireaclib, icoeff

    allocate( ctemp_rate(7, number_reaclib_sets) )
    allocate( rate_start_idx(nrat_reaclib) )
    allocate( rate_extra_mult(nrat_reaclib) )

    open(newunit=unit, file='reaclib_rate_metadata.dat')

    do ireaclib = 1, number_reaclib_sets
       do icoeff = 1, 7
          read(unit, *) ctemp_rate(icoeff, ireaclib)
       enddo
    enddo

    do ireaclib = 1, nrat_reaclib
       read(unit, *) rate_start_idx(ireaclib)
    enddo

    do ireaclib = 1, nrat_reaclib
       read(unit, *) rate_extra_mult(ireaclib)
    enddo

    close(unit)

    !$acc update device(ctemp_rate, rate_start_idx, rate_extra_mult)

  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate )
    deallocate( rate_start_idx )
    deallocate( rate_extra_mult )
  end subroutine term_reaclib


  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))


    call screening_init()
  end subroutine net_screening_init


  subroutine net_screening_finalize()
    ! Call screening_finalize

    call screening_finalize()

  end subroutine net_screening_finalize


  subroutine reaclib_evaluate(pstate, temp, iwhich, reactvec)
    !$acc routine seq

    implicit none

    type(plasma_state), intent(in) :: pstate
    real(rt), intent(in) :: temp
    integer, intent(in) :: iwhich

    real(rt), intent(inout) :: reactvec(num_rate_groups+2)
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

    real(rt) :: rate, scor ! Rate and Screening Factor
    real(rt) :: drate_dt, dscor_dt ! Temperature derivatives
    real(rt) :: dscor_dd
    real(rt) :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m, istart

    !$gpu

    ri = 0.0d0
    rate = 0.0d0
    drate_dt = 0.0d0
    irate = 0.0d0
    dirate_dt = 0.0d0
    T9 = temp/1.0d9
    T9_exp = 0.0d0

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

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = 1.0d0
    reactvec(i_dscor_dt) = 0.0d0
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

end module reaclib_rates
