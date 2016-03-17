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
    YDOT(net_meta%ip) = ( &
       - dens * Y(net_meta%ip) * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
       - dens * Y(net_meta%ic13) * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
       - dens * Y(net_meta%ip) * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
       - dens * Y(net_meta%ip) * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
       - dens * Y(net_meta%ip) * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
       )

    YDOT(net_meta%ihe4) = ( &
       + dens * Y(net_meta%ip) * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
       )

    YDOT(net_meta%ic12) = ( &
       - dens * Y(net_meta%ip) * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
       + dens * Y(net_meta%ip) * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
       )

    YDOT(net_meta%ic13) = ( &
       - dens * Y(net_meta%ic13) * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
       + Y(net_meta%in13) * rxn_rates(net_meta%k_n13_c13) &
       )

    YDOT(net_meta%in13) = ( &
       - Y(net_meta%in13) * rxn_rates(net_meta%k_n13_c13) &
       - dens * Y(net_meta%ip) * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
       + dens * Y(net_meta%ip) * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
       )

    YDOT(net_meta%in14) = ( &
       - dens * Y(net_meta%ip) * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
       + dens * Y(net_meta%ic13) * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
       + Y(net_meta%io14) * rxn_rates(net_meta%k_o14_n14) &
       )

    YDOT(net_meta%in15) = ( &
       - dens * Y(net_meta%ip) * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
       + Y(net_meta%io15) * rxn_rates(net_meta%k_o15_n15) &
       )

    YDOT(net_meta%io14) = ( &
       - Y(net_meta%io14) * rxn_rates(net_meta%k_o14_n14) &
       + dens * Y(net_meta%ip) * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
       )

    YDOT(net_meta%io15) = ( &
       - Y(net_meta%io15) * rxn_rates(net_meta%k_o15_n15) &
       + dens * Y(net_meta%ip) * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
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

      DJAC(net_meta%ip,net_meta%ip) = ( &
         - dens * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
         - dens * Y(net_meta%ic13) * rxn_rates(net_meta%k_c13_pg_n14) &
         - dens * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
         - dens * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
         - dens * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
         )

      DJAC(net_meta%ip,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ic12) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c12_pg_n13) &
         )

      DJAC(net_meta%ip,net_meta%ic13) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
         )

      DJAC(net_meta%ip,net_meta%in13) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n13_pg_o14) &
         )

      DJAC(net_meta%ip,net_meta%in14) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n14_pg_o15) &
         )

      DJAC(net_meta%ip,net_meta%in15) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%ip,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ip) = ( &
         + dens * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%ihe4,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%in13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%in15) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%ihe4,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ihe4,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ip) = ( &
         - dens * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
         + dens * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%ic12,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%ic12) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c12_pg_n13) &
         )

      DJAC(net_meta%ic12,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%in13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%in15) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%ic12,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic12,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%ip) = ( &
         - dens * Y(net_meta%ic13) * rxn_rates(net_meta%k_c13_pg_n14) &
         )

      DJAC(net_meta%ic13,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%ic13) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
         )

      DJAC(net_meta%ic13,net_meta%in13) = ( &
         +    rxn_rates(net_meta%k_n13_c13) &
         )

      DJAC(net_meta%ic13,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%in15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ic13,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%ip) = ( &
         - dens * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
         + dens * Y(net_meta%ic12) * rxn_rates(net_meta%k_c12_pg_n13) &
         )

      DJAC(net_meta%in13,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%ic12) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c12_pg_n13) &
         )

      DJAC(net_meta%in13,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%in13) = ( &
         -    rxn_rates(net_meta%k_n13_c13) &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n13_pg_o14) &
         )

      DJAC(net_meta%in13,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%in15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in13,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in14,net_meta%ip) = ( &
         - dens * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
         + dens * Y(net_meta%ic13) * rxn_rates(net_meta%k_c13_pg_n14) &
         )

      DJAC(net_meta%in14,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in14,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in14,net_meta%ic13) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_c13_pg_n14) &
         )

      DJAC(net_meta%in14,net_meta%in13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in14,net_meta%in14) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n14_pg_o15) &
         )

      DJAC(net_meta%in14,net_meta%in15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in14,net_meta%io14) = ( &
         +    rxn_rates(net_meta%k_o14_n14) &
         )

      DJAC(net_meta%in14,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%ip) = ( &
         - dens * Y(net_meta%in15) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%in15,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%in13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%in15) = ( &
         - dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n15_pa_c12) &
         )

      DJAC(net_meta%in15,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in15,net_meta%io15) = ( &
         +    rxn_rates(net_meta%k_o15_n15) &
         )

      DJAC(net_meta%io14,net_meta%ip) = ( &
         + dens * Y(net_meta%in13) * rxn_rates(net_meta%k_n13_pg_o14) &
         )

      DJAC(net_meta%io14,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io14,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io14,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io14,net_meta%in13) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n13_pg_o14) &
         )

      DJAC(net_meta%io14,net_meta%in14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io14,net_meta%in15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io14,net_meta%io14) = ( &
         -    rxn_rates(net_meta%k_o14_n14) &
         )

      DJAC(net_meta%io14,net_meta%io15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%ip) = ( &
         + dens * Y(net_meta%in14) * rxn_rates(net_meta%k_n14_pg_o15) &
         )

      DJAC(net_meta%io15,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%ic12) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%ic13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%in13) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%in14) = ( &
         + dens * Y(net_meta%ip) * rxn_rates(net_meta%k_n14_pg_o15) &
         )

      DJAC(net_meta%io15,net_meta%in15) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%io14) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%io15,net_meta%io15) = ( &
         -    rxn_rates(net_meta%k_o15_n15) &
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
    G(1) = Y(net_meta%ip)
    IER = 0
  end subroutine FCVROOTFN

end module network
