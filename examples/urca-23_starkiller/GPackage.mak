F90sources += actual_network.F90
F90sources += reaclib_rates.F90
F90sources += table_rates.F90
f90sources += physical_constants.f90

ifneq ($(USE_REACT), FALSE)
  F90sources += actual_burner.F90
  F90sources += actual_rhs.F90

  USE_SCREENING = TRUE
  USE_NEUTRINOS = TRUE
endif
