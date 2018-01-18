program call_table
  use table_rates

  implicit none

  double precision :: rhoy, temp, fret
  integer :: which_table, which_var
  logical :: continue_input = .true.

  call init_table_meta()
  
  do while ( continue_input )
     write(*,*) 'Query the tables as: <table number> <logtemp value> <logrhoy value> <table index>'
     read(*,*) which_table, temp, rhoy, which_var
     call table_meta(which_table)%bl_lookup(10.0**rhoy, 10.0**temp, which_var, fret)
     write(*,*) fret
     write(*,*) 'To continue, enter 1, to exit, enter anything else.'
     read(*,*) which_table
     if ( .not. which_table .eq. 1 ) then
        continue_input = .false.
     end if
  end do

  call term_table_meta()
end program call_table
