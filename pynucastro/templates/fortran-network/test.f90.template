program test

  use pynet

  implicit none

  double precision :: rho, T
  double precision :: X(nspec), dYdt(nspec)
  double precision :: Jac(nspec, nspec)
  double precision :: enuc

  integer :: i, j

  rho = 2.e8
  T = 1.e9

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
