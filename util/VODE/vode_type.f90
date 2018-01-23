module vode_type_module

  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: VODE_NEQS = neqs

contains

  subroutine clean_state(y, rpar)

    use bl_types, only: dp_t
    use bl_constants_module, only: ONE
    use network, only: nspec_evolve
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps
    use extern_probin_module, only: renormalize_abundances, small_x_safe

    implicit none

    real(dp_t) :: y(neqs), rpar(n_rpar_comps)

    real(dp_t) :: small_temp

    ! Ensure that mass fractions always stay positive and less than or equal to 1.

    y(1:nspec_evolve) = max(min(y(1:nspec_evolve), ONE), small_x_safe)

    ! Renormalize the abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

  end subroutine clean_state


  subroutine renormalize_species(y, rpar)

    use bl_types, only: dp_t
    use network, only: aion, aion_inv, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps, irp_nspec, n_not_evolved

    implicit none

    real(dp_t) :: y(neqs), rpar(n_rpar_comps)

    real(dp_t) :: nspec_sum

    nspec_sum = &
         sum(y(1:nspec_evolve)) + &
         sum(rpar(irp_nspec:irp_nspec+n_not_evolved-1))

    y(1:nspec_evolve) = y(1:nspec_evolve) / nspec_sum
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1) / nspec_sum

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine burn_to_vode(state, y, rpar, ydot, jac)

    use bl_types, only: dp_t
    use network, only: nspec, nspec_evolve
    use rpar_indices, only: irp_dens, irp_temp, irp_nspec, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, burn_t, net_ienuc

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(neqs)
    real(dp_t), optional :: ydot(neqs), jac(neqs, neqs)

    integer :: n

    rpar(irp_dens) = state % rho
    rpar(irp_temp) = state % T

    y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = state % xn(nspec_evolve+1:nspec)

    y(net_ienuc) = state % e

    if (present(ydot)) then
       ydot = state % ydot
    endif

    if (present(jac)) then
       jac = state % jac
    endif

  end subroutine burn_to_vode



  ! Given an rpar array and the integration state, set up a burn state.

  subroutine vode_to_burn(y, rpar, state)

    use bl_types, only: dp_t
    use network, only: nspec, nspec_evolve
    use rpar_indices, only: irp_dens, irp_temp, irp_nspec, n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, burn_t, net_ienuc

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(neqs)

    integer :: n

    state % rho      = rpar(irp_dens)
    state % T        = rpar(irp_temp)
    state % e        = y(net_ienuc)

    state % xn(1:nspec_evolve) = y(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1)

  end subroutine vode_to_burn

end module vode_type_module
