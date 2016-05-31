f90sources += actual_network.f90
f90sources += actual_network_data.f90
f90sources += burn_type.f90
f90sources += net_rates.f90
f90sources += table_rates.f90
f90sources += physical_constants.f90
f90sources += terminator.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_burner_data.f90
  f90sources += actual_rhs.f90
  f90sources += sneut5.f90
  f90sources += terminator.f90

  USE_SCREENING = TRUE
endif
