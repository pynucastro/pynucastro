! This program writes back the tabulated rates for
! reproducing the plots of, e.g. Toki, et al 2015.

program output_table
  use table_rates

  implicit none

  double precision :: lrhoy, ltemp, lrhoy_step, ltemp_step, fret
  integer :: ilrhoy, iltemp
  integer :: which_table, which_var
  character(len=50) :: rfmt
  character(len=50) :: hfmt
  double precision, parameter :: lrhoy_start = 8.0d0
  double precision, parameter :: lrhoy_end   = 9.2d0
  double precision, parameter :: ltemp_start = 8.0d0
  double precision, parameter :: ltemp_end   = 9.2d0 
  integer, parameter :: lrhoy_num = 1000
  integer, parameter :: ltemp_num = 7
  integer, parameter :: NPROFILE = 3

  lrhoy_step = (lrhoy_end - lrhoy_start)/(lrhoy_num - 1)
  ltemp_step = (ltemp_end - ltemp_start)/(ltemp_num - 1)
  
  write(hfmt,'(A,I5,A)') '(', NPROFILE, '(A18))'
  write(rfmt,'(A,I5,A)') '(', NPROFILE, '(ES25.14))'

  which_var = jtab_rate

  ! Initialize tables
  call init_table_meta()

  ! Emission
  which_table = 1
  open(unit=1, file='output_emission_table.dat', recl=(25*NPROFILE+10), form='formatted')
  write(1,*) ltemp_num
  write(1,*) lrhoy_num
  write(1, fmt=hfmt) 'log_lrhoy', 'log_ltemp', 'rate'
  do iltemp = 0, ltemp_num-1
     ltemp = ltemp_start + dble(iltemp)*ltemp_step
     write(*,*) ltemp
     do ilrhoy = 0, lrhoy_num-1
        lrhoy = lrhoy_start + dble(ilrhoy)*lrhoy_step
        call table_meta(which_table)%bl_lookup(10.0**lrhoy, 10.0**ltemp, which_var, fret)
        write(1, fmt=rfmt) lrhoy, ltemp, fret
     end do
  end do
  close(unit=1)

  ! Capture
  which_table = 2
  open(unit=2, file='output_capture_table.dat', recl=(25*NPROFILE+10), form='formatted')
  write(2,*) ltemp_num
  write(2,*) lrhoy_num
  write(2, fmt=hfmt) 'log_lrhoy', 'log_ltemp', 'rate'
  do iltemp = 0, ltemp_num-1
     ltemp = ltemp_start + dble(iltemp)*ltemp_step
     write(*,*) ltemp
     do ilrhoy = 0, lrhoy_num-1
        lrhoy = lrhoy_start + dble(ilrhoy)*lrhoy_step
        call table_meta(which_table)%bl_lookup(10.0**lrhoy, 10.0**ltemp, which_var, fret)
        write(2, fmt=rfmt) lrhoy, ltemp, fret
     end do
  end do
  close(unit=2)

  ! Deallocate table memory
  call term_table_meta()
end program output_table
