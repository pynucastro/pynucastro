module burn_type_module

  use bl_types, only: dp_t
  use network, only: nspec, nspec_evolve

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn.

  ! Set the number of independent variables -- this should be
  ! temperature, enuc + the number of species which participate
  ! in the evolution equations.

  integer, parameter :: neqs = 1 + nspec_evolve

  ! Indices of the temperature and energy variables in the work arrays.

  integer, parameter :: net_ienuc = nspec_evolve + 1

  type :: burn_t

    real(dp_t) :: rho
    real(dp_t) :: T
    real(dp_t) :: e
    real(dp_t) :: xn(nspec)

    real(dp_t) :: ye

    ! The following are the actual integration data.
    ! To avoid potential incompatibilities we won't
    ! include the integration vector y itself here.
    ! It can be reconstructed from all of the above
    ! data, particularly xn, e, and T.

    real(dp_t) :: ydot(neqs)
    real(dp_t) :: jac(neqs, neqs)

    ! diagnostics
    integer :: n_rhs
    integer :: n_jac

    ! Integration time.

    real(dp_t) :: time

  end type burn_t

contains

  subroutine normalize_abundances_burn(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use probin_module, only: small_x

    implicit none

    type (burn_t), intent(inout) :: state

    state % xn(:) = max(small_x, min(ONE, state % xn(:)))
    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances_burn

  subroutine update_electron_fraction(state)

    !$acc routine seq

    use network, only: zion, aion_inv
    
    implicit none

    type (burn_t), intent(inout) :: state

    call normalize_abundances_burn(state)
    state % ye = sum(state % xn(:)*zion(:)*aion_inv(:))
    
  end subroutine update_electron_fraction

end module burn_type_module
