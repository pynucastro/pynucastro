module pynet
    integer, parameter :: nspec = 9

    double precision, parameter :: aion(nspec) = [ &
        1, 1, 4, 12,  &
        16, 20, 23, 23,  &
        24]

    double precision, parameter :: zion(nspec) = [ &
        0, 1, 2, 6,  &
        8, 10, 11, 12,  &
        12]

    character(len=6), dimension(nspec) :: spec_names = [ &
        character(len=6) :: "n", "h1", "he4", "c12",  &
                            "o16", "ne20", "na23", "mg23",  &
                            "mg24"]

    interface
       subroutine rhs_f(rho, T, X, dYdt, size) bind(C, name="rhs_f")
         use, intrinsic :: iso_c_binding, only: c_double, c_int
         import nspec
         real(kind=c_double), intent(in), value :: rho, T
         real(kind=c_double), intent(in), dimension(nspec) :: X
         real(kind=c_double), intent(inout), dimension(nspec) :: dYdt
         integer(c_int), value :: size
       end subroutine rhs_f

       subroutine jac_f(rho, T, X, J, size) bind(C, name="jac_f")
         use, intrinsic :: iso_c_binding, only: c_double, c_int
         import nspec
         real(kind=c_double), intent(in), value :: rho, T
         real(kind=c_double), intent(in), dimension(nspec) :: X
         real(kind=c_double), intent(inout), dimension(nspec, nspec) :: J
         integer(c_int), value :: size
       end subroutine jac_f

       subroutine ener_gener_f(dYdt, size, enuc) bind(C, name="ener_gener_f")
         use, intrinsic :: iso_c_binding, only: c_double, c_int
         import nspec
         real(kind=c_double), intent(in), dimension(nspec) :: dYdt
         integer(c_int), value :: size
         real(kind=c_double), intent(out) :: enuc
       end subroutine ener_gener_f

       subroutine network_init() bind(C, name="network_init")
       end subroutine network_init

    end interface

end module pynet


