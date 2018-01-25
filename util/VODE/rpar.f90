! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module rpar_indices

  use network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_temp = irp_dens + 1
  integer, parameter :: irp_nspec = irp_temp + 1
  integer, parameter :: irp_y_init = irp_nspec + n_not_evolved
  integer, parameter :: irp_t0 = irp_y_init + neqs
  integer, parameter :: n_rpar_comps = irp_t0

end module rpar_indices
