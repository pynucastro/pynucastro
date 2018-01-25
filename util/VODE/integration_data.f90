module integration_data

  use bl_types, only: dp_t

  implicit none

  type :: integration_status_t

     logical :: integration_complete
     real(dp_t) :: atol_spec, atol_enuc
     real(dp_t) :: rtol_spec, rtol_enuc

  end type integration_status_t

end module integration_data
