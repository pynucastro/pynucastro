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

    integer, parameter :: n_ = 1
    integer, parameter :: h1_ = 2
    integer, parameter :: he4_ = 3
    integer, parameter :: c12_ = 4
    integer, parameter :: o16_ = 5
    integer, parameter :: ne20_ = 6
    integer, parameter :: na23_ = 7
    integer, parameter :: mg23_ = 8
    integer, parameter :: mg24_ = 9

    interface
       subroutine rhs_f(rho, T, X, dYdt) bind(C, name="rhs_f")
         use, intrinsic :: iso_c_binding, only: c_double
         import nspec
         real(kind=c_double), intent(in), value :: rho, T
         real(kind=c_double), intent(in), dimension(nspec) :: X
         real(kind=c_double), intent(inout), dimension(nspec) :: dYdt
       end subroutine rhs_f

       subroutine jac_f(rho, T, X, J) bind(C, name="jac_f")
         use, intrinsic :: iso_c_binding, only: c_double
         import nspec
         real(kind=c_double), intent(in), value :: rho, T
         real(kind=c_double), intent(in), dimension(nspec) :: X
         real(kind=c_double), intent(inout), dimension(nspec, nspec) :: J
       end subroutine jac_f

       subroutine ener_gener_f(dYdt, enuc) bind(C, name="ener_gener_f")
         use, intrinsic :: iso_c_binding, only: c_double
         import nspec
         real(kind=c_double), intent(in), dimension(nspec) :: dYdt
         real(kind=c_double), intent(out) :: enuc
       end subroutine ener_gener_f

       subroutine network_init() bind(C, name="network_init")
       end subroutine network_init

    end interface

end module pynet


