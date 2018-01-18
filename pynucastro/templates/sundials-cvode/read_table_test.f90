! This program tests table-reading and writing
program read_table_test

  implicit none

  double precision, target, dimension(:,:,:), allocatable :: rate_table_temp
  double precision, target, dimension(:,:,:), allocatable :: rate_table_1
  character(len=*), parameter :: rate_table_file_1 = '23Ne-23Na_eemission.dat'
  character(len=*), parameter :: rate_table_test_1 = '23Ne-23Na_eemission.out.min'
  integer, parameter :: num_header = 7
  integer, parameter :: num_rhoy = 152
  integer, parameter :: num_temp = 39
  integer, parameter :: num_vars = 6
  integer :: i, j, k

  allocate( rate_table_temp(num_temp, num_rhoy, num_vars+2) )
  
  open(unit=11, file=rate_table_file_1)
  do i = 1, num_header
     read(11,*)
  end do
  do j = 1, num_rhoy
     do i = 1, num_temp
        read(11,*) ( rate_table_temp(i, j, k), k=1, num_vars+2 )
     end do
     if (j/=num_rhoy) then
        read(11,*)
     end if
  end do
  close(11)

  allocate ( rate_table_1(num_temp, num_rhoy, num_vars) )
  rate_table_1(:,:,:) = rate_table_temp(:,:,3:num_vars+2)

  open(unit=12, file=rate_table_test_1)
  do j = 1, num_rhoy
     do i = 1, num_temp
        write(12, fmt='(8 (ES25.14))') ( rate_table_1(i, j, k), k=1, num_vars )
     end do
  end do
  close(12)

end program read_table_test
