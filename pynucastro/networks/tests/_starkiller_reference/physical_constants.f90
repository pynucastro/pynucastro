module physical_constants

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter :: pi = 3.141592653589793238462643383279502884197e0_rt
  real(rt), parameter :: cm_per_fm = 1.0e-13_rt

  ! From 2014 CODATA recommended values
  real(rt), parameter :: clight  = 299792458e2_rt  ! cm/s

  real(rt), parameter :: N_avo = 6.022140857e23_rt ! 1/mol

  real(rt), parameter :: gram_per_amu = 1660539.040e-30_rt

  real(rt), parameter :: erg_per_eV  = 1.6021766208e-12_rt
  real(rt), parameter :: erg_per_MeV = erg_per_eV * 1.0e6_rt
  real(rt), parameter :: gram_per_MeV = erg_per_MeV / clight**2

end module physical_constants
