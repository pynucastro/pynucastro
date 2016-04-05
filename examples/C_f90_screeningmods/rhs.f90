module rhs
  use physical_constants, only: N_AVO
  use cvode_parameters, only: cv_pars
  use network
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
   
    
    write(*,*) "T: ", T
    YDOT(jn) = ( &
       - Y(jn) * rxn_rates(3, k_n_p) * rxn_rates(1, k_n_p) &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
       )

    YDOT(jp) = ( &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
       + Y(jn) * rxn_rates(3, k_n_p) * rxn_rates(1, k_n_p) &
       )

    YDOT(jhe4) = ( &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
       )

    YDOT(jc12) = ( &
       - 2 * 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
       - 2 * 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
       - 2 * 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
       )

    YDOT(jne20) = ( &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
       )

    YDOT(jna23) = ( &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
       )

    YDOT(jmg23) = ( &
       + 0.500000000000000 * dens * Y(jc12)**2 * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
       )


    YDOT(jenuc) = 0.0d0
    do i = 1, nspec
      YDOT(jenuc) = YDOT(jenuc) + N_AVO * ebind(i) * YDOT(i)
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

       if ( screen_reaclib ) then
         call fill_plasma_state(state, temp, dens, Y(1:nspec))
       end if
    
       do i = 1, nreact
         call rate_evaluate(state, rhoy, cv_pars%temp, i, rxn_rates(:,i))
       end do

      ! DJAC(j, i) = d(YDOT(j))/dY(i)

      DJAC(jn,jn) = ( &
         -    rxn_rates(3, k_n_p) * rxn_rates(1, k_n_p) &
         )

      DJAC(jn,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jn,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jn,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
         )

      DJAC(jn,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jn,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jn,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jp,jn) = ( &
         +    rxn_rates(3, k_n_p) * rxn_rates(1, k_n_p) &
         )

      DJAC(jp,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jp,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jp,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
         )

      DJAC(jp,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jp,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jp,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jn) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
         )

      DJAC(jhe4,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jhe4,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jn) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jc12) = ( &
         - 2 * 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
         - 2 * 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
         - 2 * 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
         )

      DJAC(jc12,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jc12,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jn) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12a_ne20) * rxn_rates(1, k_c12_c12a_ne20) &
         )

      DJAC(jne20,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jne20,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jn) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12p_na23) * rxn_rates(1, k_c12_c12p_na23) &
         )

      DJAC(jna23,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jna23,jmg23) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jn) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jp) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jhe4) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jc12) = ( &
         + 0.500000000000000 * dens * 2*Y(jc12) * rxn_rates(3, k_c12_c12n_mg23) * rxn_rates(1, k_c12_c12n_mg23) &
         )

      DJAC(jmg23,jne20) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jna23) = ( &
         + 0.0d0 &
         )

      DJAC(jmg23,jmg23) = ( &
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
