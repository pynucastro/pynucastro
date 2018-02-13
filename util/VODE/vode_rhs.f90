  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use network, only: aion, nspec_evolve
    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc
    use bl_constants_module, only: ZERO, ONE
    use rhs_module, only: netrhs
    use vode_type_module, only: clean_state, renormalize_species, &
                                burn_to_vode, vode_to_burn
    use rpar_indices, only: n_rpar_comps, irp_y_init

    implicit none

    integer,    intent(IN   ) :: neq, ipar
    real(dp_t), intent(INOUT) :: time, y(neq)
    real(dp_t), intent(INOUT) :: rpar(n_rpar_comps)
    real(dp_t), intent(  OUT) :: ydot(neq)

    type (burn_t) :: burn_state

    real(dp_t) :: limit_factor, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_ienuc) = denuc/dt

    ydot = ZERO

    ! Fix the state as necessary.

    call clean_state(y, rpar)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(y, rpar, burn_state)

    burn_state % time = time
    call netrhs(burn_state)

    ! We integrate X, not Y
    burn_state % ydot(1:nspec_evolve) = &
         burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    call burn_to_vode(burn_state, y, rpar, ydot = ydot)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    use network, only: aion, aion_inv, nspec_evolve
    use bl_constants_module, only: ZERO, ONE
    use rhs_module, only: netjac
    use burn_type_module, only: burn_t, net_ienuc
    use vode_type_module, only: vode_to_burn, burn_to_vode
    use rpar_indices, only: n_rpar_comps, irp_y_init
    use bl_types, only: dp_t

    implicit none

    integer   , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    real(dp_t), intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
    real(dp_t), intent(  OUT) :: pd(neq,neq)

    type (burn_t) :: state
    real(dp_t) :: limit_factor, t_enuc
    integer :: n

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    state % time = time
    call netjac(state)

    ! We integrate X, not Y
    do n = 1, nspec_evolve
       state % jac(n,:) = state % jac(n,:) * aion(n)
       state % jac(:,n) = state % jac(:,n) * aion_inv(n)
    enddo

    call burn_to_vode(state, y, rpar, jac = pd)

  end subroutine jac
