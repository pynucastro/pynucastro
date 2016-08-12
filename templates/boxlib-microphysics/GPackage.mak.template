f90sources += actual_network.f90
f90sources += reaclib_rates.f90
f90sources += table_rates.f90
f90sources += physical_constants.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_rhs.f90

  USE_SCREENING = TRUE
  USE_NEUTRINOS = TRUE
endif
