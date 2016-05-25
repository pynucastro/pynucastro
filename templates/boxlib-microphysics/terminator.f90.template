module terminator_module

  implicit none

contains
  
  subroutine terminate()
    use net_rates, only: term_reaclib
    use table_rates, only: term_table_meta

    implicit none
    
    call term_reaclib()
    call term_table_meta()
  end subroutine terminate
  
end module terminator_module
