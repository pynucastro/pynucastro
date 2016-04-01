module rhs
  use physical_constants, only: pc
  use cvode_parameters, only: cv_pars
  use network, only: nspec, ebind, net_meta
  use net_rates, only: nreact, screen_reaclib, rate_evaluate
  use screening_module, only: plasma_state, fill_plasma_state

  implicit none

contains

  subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER) bind(C, name='fcvfun_')
    type(plasma_state) :: state
    double precision, dimension(*), intent(inout) :: Y, RPAR
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(out) :: YDOT
    double precision, intent(in) :: T
    integer, intent(out) :: IER

    double precision, dimension(4, nreact) :: rxn_rates
    integer :: i
    double precision :: dens, temp, ye, rhoy

    IER = -1

    ! Enforce all Y's are positive
    do i = 1, nspec
       Y(i) = max(0.0d0, Y(i))
    end do
    
    ! Calculate rates
    ! density is cv_pars%dens
    ! temperature is cv_pars%temp

    dens = cv_pars%dens
    temp = cv_pars%temp
    ye   = cv_pars%ye
    rhoy = dens*ye

    call fill_plasma_state(state, temp, dens, Y(1:nspec))
    do i = 1, nreact
      call rate_evaluate(state, rhoy, cv_pars%temp, i, rxn_rates(:,i))
    end do
   
    ! Add screening factors to the rates
    
    write(*,*) "T: ", T
    YDOT(net_meta%in) = ( &
       - Y(net_meta%in) * rxn_rates(3,net_meta%k_n_p) * rxn_rates(1,net_meta%k_n_p) &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
       )

    YDOT(net_meta%ip) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
       + Y(net_meta%in) * rxn_rates(3,net_meta%k_n_p) * rxn_rates(1,net_meta%k_n_p) &
       )

    YDOT(net_meta%ihe4) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
       )

    YDOT(net_meta%ic12) = ( &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
       - 2 * 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
       )

    YDOT(net_meta%ine20) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
       )

    YDOT(net_meta%ina23) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
       )

    YDOT(net_meta%img23) = ( &
       + 0.500000000000000 * dens * Y(net_meta%ic12)**2 * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
       )


    YDOT(net_meta%ienuc) = 0.0d0
    do i = 1, nspec
      YDOT(net_meta%ienuc) = YDOT(net_meta%ienuc) + pc%N_Avogadro * ebind(i) * YDOT(i)
    end do

    write(*,*) '______________________________'
    do i = 1, nspec
      write(*,*) 'YDOT(',i,'): ',YDOT(i)
    end do

    IER = 0 ! Successful
  end subroutine FCVFUN

  subroutine FCVDJAC(NEQ, T, Y, FY, DJAC, H_STEP, IPAR, RPAR, WK1, WK2, WK3, IER) bind(C, name='fcvdjac_')
    type(plasma_state) :: state
    integer, intent(in) :: NEQ ! number of ODEs
    double precision, intent(in) :: T ! independent variable
    double precision, dimension(*), intent(inout) :: Y, FY ! y and its derivative
    double precision, dimension(NEQ,*), intent(out) :: DJAC ! dense Jacobian
    double precision, intent(in) :: H_STEP ! current stepsize
    integer, dimension(*), intent(in) :: IPAR
    double precision, dimension(*), intent(in) :: RPAR
    double precision, dimension(NEQ), intent(in) :: WK1, WK2, WK3
    integer, intent(out) :: IER

    double precision, dimension(4, nreact) :: rxn_rates
    double precision :: dens, temp, ye, rhoy
    integer :: i, j

    IER = -1

    ! Enforce all Y's are positive
    do i = 1, nspec
       Y(i) = max(0.0d0, Y(i))
    end do

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
       ! Add screening to the jacobian
       if ( screen_reaclib ) then
         call fill_plasma_state(state, temp, dens, Y(1:nspec))
       end if
    
       do i = 1, nreact
         call rate_evaluate(state, rhoy, cv_pars%temp, i, rxn_rates(:,i))
       end do

      ! DJAC(j, i) = d(YDOT(j))/dY(i)

      DJAC(net_meta%in,net_meta%in) = ( &
         -    rxn_rates(3,net_meta%k_n_p) * rxn_rates(1,net_meta%k_n_p) &
         )

      DJAC(net_meta%in,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
         )

      DJAC(net_meta%in,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%in,net_meta%img23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%in) = ( &
         +    rxn_rates(3,net_meta%k_n_p) * rxn_rates(1,net_meta%k_n_p) &
         )

      DJAC(net_meta%ip,net_meta%ip) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ihe4) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ip,net_meta%ic12) = ( &
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ip,net_meta%ine20) = ( &
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
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
         )

      DJAC(net_meta%ihe4,net_meta%ine20) = ( &
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
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
         - 2 * 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ic12,net_meta%ine20) = ( &
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
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12a_ne20) * rxn_rates(1,net_meta%k_c12_c12a_ne20) &
         )

      DJAC(net_meta%ine20,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%ina23) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ine20,net_meta%img23) = ( &
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
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12p_na23) * rxn_rates(1,net_meta%k_c12_c12p_na23) &
         )

      DJAC(net_meta%ina23,net_meta%ine20) = ( &
         + 0.0d0 &
         )

      DJAC(net_meta%ina23,net_meta%ina23) = ( &
         + 0.0d0 &
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
         + 0.500000000000000 * dens * 2*Y(net_meta%ic12) * rxn_rates(3,net_meta%k_c12_c12n_mg23) * rxn_rates(1,net_meta%k_c12_c12n_mg23) &
         )

      DJAC(net_meta%img23,net_meta%ine20) = ( &
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

end module rhs
