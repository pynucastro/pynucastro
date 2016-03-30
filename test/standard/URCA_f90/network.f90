module network

  use physical_constants
  use cvode_parameters
  use net_rates

  implicit none

contains

  subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER) bind(C, name='fcvfun_')
    double precision, dimension(*), intent(in) :: Y, RPAR
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(out) :: YDOT
    double precision, intent(in) :: T
    integer, intent(out) :: IER

    double precision, dimension(net_meta%nrxn) :: rxn_rates
    integer :: i
    double precision :: ri, dens, temp, ye, rhoy

    IER = -1

    ! Calculate rates
    ! density is cv_pars%dens
    ! temperature is cv_pars%temp

    dens = cv_pars%dens
    temp = cv_pars%temp
    ye   = cv_pars%ye
    rhoy = dens*ye

    ri = 0.0d0
    do i = 1, net_meta%nrxn
      call net_meta%evaluate(rhoy, cv_pars%temp, i, ri)
      rxn_rates(i) = ri
    end do

    write(*,*) "T: ", T
    YDOT(net_meta%in) = ( &
       - Y(net_meta%in) * rxn_rates(net_meta%k_n_p) &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12n_mg23) &
       )

    YDOT(net_meta%ip) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12p_na23) &
       + Y(net_meta%in) * rxn_rates(net_meta%k_n_p) &
       )

    YDOT(net_meta%ihe4) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12a_ne20) &
       )

    YDOT(net_meta%ic12) = ( &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12a_ne20) &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12n_mg23) &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12p_na23) &
       )

    YDOT(net_meta%ine20) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12a_ne20) &
       )

    YDOT(net_meta%ine23) = ( &
       - Y(net_meta%ine23) * rxn_rates(net_meta%k_ne23_na23) &
       + Y(net_meta%ina23) * rxn_rates(net_meta%k_na23_ne23) &
       )

    YDOT(net_meta%ina23) = ( &
       - Y(net_meta%ina23) * rxn_rates(net_meta%k_na23_ne23) &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12p_na23) &
       + Y(net_meta%ine23) * rxn_rates(net_meta%k_ne23_na23) &
       )

    YDOT(net_meta%img23) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(net_meta%k_c12_c12n_mg23) &
       )


    YDOT(net_meta%ienuc) = 0.0d0
    do i = 1, net_meta%nnuc
      YDOT(net_meta%ienuc) = YDOT(net_meta%ienuc) + pc%N_Avogadro * net_meta%ebind(i) * YDOT(i)
    end do

    write(*,*) '______________________________'
    do i = 1, net_meta%nnuc
      write(*,*) 'YDOT(',i,'): ',YDOT(i)
    end do

    IER = 0 ! Successful
  end subroutine FCVFUN

  subroutine FCVDJAC(NEQ, T, Y, FY, DJAC, H_STEP, IPAR, RPAR, WK1, WK2, WK3, IER) bind(C, name='fcvdjac_')
    integer, intent(in) :: NEQ ! number of ODEs
    double precision, intent(in) :: T ! independent variable
    double precision, dimension(*), intent(in) :: Y, FY ! y and its derivative
    double precision, dimension(NEQ,*), intent(out) :: DJAC ! dense Jacobian
    double precision, intent(in) :: H_STEP ! current stepsize
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    double precision, dimension(NEQ), intent(in) :: WK1, WK2, WK3
    integer, intent(out) :: IER

    double precision, dimension(net_meta%nrxn) :: rxn_rates
    double precision :: ri, dens, temp, ye, rhoy
    integer :: i, j

    IER = -1

    dens = cv_pars%dens
    temp = cv_pars%temp
    ye   = cv_pars%ye
    rhoy = dens*ye

    if (T == 0.0d0) then
    
      do i = 1, NEQ
        do j = 1, NEQ
          DJAC(j, i) = 0.0d0
        end do
      end do

    else
    
       ri = 0.0d0

       do i = 1, net_meta%nrxn
          call net_meta%evaluate(rhoy, temp, i, ri)
          rxn_rates(i) = ri
       end do

      ! DJAC(j, i) = d(YDOT(j))/dY(i)

      DJAC(net_meta%in,net_meta%in) = ( &
         -    rxn_rates(net_meta%k_n_p) &
         )

      DJAC(net_meta%in,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12n_mg23) &
         )

      DJAC(net_meta%in,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%in) = ( &
         +    rxn_rates(net_meta%k_n_p) &
         )

      DJAC(net_meta%ip,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ip,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12a_ne20) &
         )

      DJAC(net_meta%ihe4,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ic12) = ( &
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12a_ne20) &
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12n_mg23) &
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ic12,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12a_ne20) &
         )

      DJAC(net_meta%ine20,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine23,net_meta%ine23) = ( &
         -    rxn_rates(net_meta%k_ne23_na23) &
         )

      DJAC(net_meta%ine23,net_meta%ina23) = ( &
         +    rxn_rates(net_meta%k_na23_ne23) &
         )

      DJAC(net_meta%ine23,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ina23,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%ine23) = ( &
         +    rxn_rates(net_meta%k_ne23_na23) &
         )

      DJAC(net_meta%ina23,net_meta%ina23) = ( &
         -    rxn_rates(net_meta%k_na23_ne23) &
         )

      DJAC(net_meta%ina23,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%in) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_c12n_mg23) &
         )

      DJAC(net_meta%img23,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%ine23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%img23,net_meta%img23) = ( &
         + 0.0d0 &
         )


    end if

    IER = 0 ! Success
  end subroutine FCVDJAC

  subroutine FCVROOTFN(T, Y, G, IPAR, RPAR, IER) bind(C, name='fcvrootfn_')
    double precision, intent(in) :: T
    double precision, dimension(*), intent(in) :: Y
    double precision, dimension(*), intent(inout) :: G
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    integer, intent(out) :: IER

    IER = -1
    ! Here you can add custom abundance stopping criteria
    G(1) = 1.0d0
    IER = 0
  end subroutine FCVROOTFN

end module network
