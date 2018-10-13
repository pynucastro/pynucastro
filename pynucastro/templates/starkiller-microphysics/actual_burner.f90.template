module actual_burner_module

  use network

  implicit none

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init
    use reaclib_rates, only: init_reaclib, net_screening_init
    use table_rates, only: init_tabular

    implicit none

    call integrator_init()
    
    call actual_network_init()

    call init_reaclib()
    call init_tabular()
    call net_screening_init()
    
  end subroutine actual_burner_init

  subroutine actual_burner_finalize
    use reaclib_rates, only: term_reaclib
    use table_rates, only: term_table_meta

    implicit none
    
    call term_reaclib()
    call term_table_meta()
  end subroutine actual_burner_finalize

  subroutine actual_burner(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator
    use burn_type_module, only: burn_t
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    real(kind=rt),  intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)
  end subroutine actual_burner

end module actual_burner_module
