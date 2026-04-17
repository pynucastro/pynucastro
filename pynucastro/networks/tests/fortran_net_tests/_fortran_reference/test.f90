program test

  use pynet

  implicit none

  double precision :: rho, T
  double precision :: X(nspec), dYdt(nspec)
  double precision :: Jac(nspec, nspec)
  double precision :: enuc

  integer :: i, j, nargs, ios
  character(len=128) :: arg

  ! we can do
  !   ./test rho T
  ! to pass in the density and temperature.
  ! composition will always be uniform.

  ! defaults
  rho = 2.e8
  T = 1.e9

  nargs = command_argument_count()
  if (nargs == 2) then
     call get_command_argument(1, arg)
     read(arg, *, iostat=ios) rho
     if (ios /= 0) stop "Error: invalid rho"

     call get_command_argument(2, arg)
     read(arg, *, iostat=ios) T
     if (ios /= 0) stop "Error: invalid T"
  else if (nargs /= 0) then
     print *, "Usage: ./test <rho> <T>"
     stop 1
  end if

  print *, "(rho, T) = ", rho, " ", T

  X(:) = 1.0 / nspec
  dYdt(:) = 0.0

  call network_init()

  call rhs_f(rho, T, X, dYdt)
  call jac_f(rho, T, X, Jac)

  do i = 1, nspec
     print *, "ydot(", spec_names(i), ") = ", dYdt(i)
  end do

  print *, " "

  do j = 1, nspec
     do i = 1, nspec
        print *, "jac(", spec_names(i), ", ", spec_names(j), ") = ", Jac(i,j)
     end do
  end do

  call ener_gener_f(dYdt, enuc)

  print *, "energy generation rate = ", enuc

end program test
