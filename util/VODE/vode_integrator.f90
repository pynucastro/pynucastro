! Common variables and routines for burners
! that use VODE for their integration.

module vode_integrator_module

  use network, only: nspec
  use rpar_indices
  use vode_type_module
  use burn_type_module
  use bl_types

  implicit none

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

  ! Tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4).
  !
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances,
  ! since we want to be easier on the temperature than the species.

  integer, parameter :: ITOL = 4

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.

  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the
  ! optional arguments flag.

  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of size 30 + NEQ. These are VODE constants
  ! that depend on the integration mode we're using -- see dvode.f.

  integer, parameter :: LRW = 22 + 9*neqs + 2*neqs**2
  integer, parameter :: LIW = 30 + neqs

contains

  subroutine vode_integrator_init()

    implicit none

  end subroutine vode_integrator_init



  ! Main interface

  subroutine vode_integrator(state_in, state_out, dt, time, status)

    use rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose
    use rhs_module, only : update_unevolved_species
    use bl_constants_module, only : ZERO, ONE
    use integration_data, only: integration_status_t

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time
    real(dp_t) :: energy_offset
    type (integration_status_t), intent(inout) :: status

    ! Local variables

    real(dp_t) :: local_time

    ! Work arrays

    real(dp_t) :: y(neqs)
    real(dp_t) :: atol(neqs), rtol(neqs)
    real(dp_t) :: rwork(LRW)
    integer    :: iwork(LIW)
    real(dp_t) :: rpar(n_rpar_comps)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar

    real(dp_t) :: sum
    real(dp_t) :: retry_change_factor

    logical :: integration_failed
    real(dp_t), parameter :: failure_tolerance = 1.d-2

    EXTERNAL jac, f_rhs

    if (jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC
    else if (jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC
    else
       write(*,*) "Error: unknown Jacobian mode in vode_integrator.f90."
       stop
    endif

    integration_failed = .false.

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec_evolve) = status % atol_spec ! mass fractions
    atol(net_ienuc)      = status % atol_enuc ! energy generated

    rtol(1:nspec_evolve) = status % rtol_spec ! mass fractions
    rtol(net_ienuc)      = status % rtol_enuc ! energy generated

    ! We want VODE to re-initialize each time we call it.

    istate = 1

    ! Initialize work arrays to zero.

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 150000

    ! Disable printing of messages about T + H == T unless we are in verbose mode.

    if (burner_verbose) then
       iwork(7) = 1
    else
       iwork(7) = 0
    endif

    ! Initialize the integration time.

    local_time = ZERO

    ! Convert our input burn state into an EOS type.

    call burn_to_vode(state_in, y, rpar)

    ! We assume that (rho, T) coming in are valid

    energy_offset = y(net_ienuc)
    y(net_ienuc) = ZERO

    ! Set the time offset -- this converts between the local integration 
    ! time and the simulation time

    rpar(irp_t0) = time

    ! Save the initial state.

    rpar(irp_y_init:irp_y_init + neqs - 1) = y

    ! Call the integration routine.

    call dvode(f_rhs, neqs, y, local_time, local_time + dt, &
               ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

    ! If we still failed, print out the current state of the integration.
    ! VODE does not always fail even though it can lead to unphysical states,
    ! so add some sanity checks that trigger a retry even if VODE thinks
    ! the integration was successful.

    if (istate < 0) then
       integration_failed = .true.
    end if

    if (any(y(1:nspec_evolve) < -failure_tolerance)) then
       integration_failed = .true.
    end if

    if (any(y(1:nspec_evolve) > 1.d0 + failure_tolerance)) then
       integration_failed = .true.
    end if

    if (integration_failed) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', local_time
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'xn current = ', y(1:nspec_evolve), rpar(irp_nspec:irp_nspec+n_not_evolved-1)
       print *, 'energy generated = ', y(net_ienuc)

       status % integration_complete = .false.
       return
    else
       status % integration_complete = .true.
    end if

    y(net_ienuc) = y(net_ienuc) + energy_offset
    
    ! Store the final data, and then normalize abundances.
    call vode_to_burn(y, rpar, state_out)

    ! Get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = iwork(12)
    state_out % n_jac = iwork(13)

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
    endif

    call normalize_abundances_burn(state_out)

    ! set the integration time for any diagnostics
    state_out % time = time + dt

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, &
                ' energy released: ', state_out % e - state_in % e
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)

    endif

  end subroutine vode_integrator

end module vode_integrator_module
