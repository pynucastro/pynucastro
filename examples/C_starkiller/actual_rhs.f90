module actual_rhs_module

  use bl_types
  use bl_constants_module
  use physical_constants, only: N_AVO
  use network
  use reaclib_rates, only: screen_reaclib, reaclib_evaluate
  use table_rates
  use screening_module, only: plasma_state, fill_plasma_state
  use sneut_module, only: sneut5
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use burn_type_module

  implicit none

  type :: rate_eval_t
     double precision :: unscreened_rates(4, nrates)
     double precision :: screened_rates(nrates)
     double precision :: dqweak(nrat_tabular)
     double precision :: epart(nrat_tabular)
  end type rate_eval_t
  
contains

  subroutine actual_rhs_init()
    ! STUB FOR MAESTRO'S TEST_REACT. ALL THE INIT IS DONE BY BURNER_INIT
    return
  end subroutine actual_rhs_init
  
  subroutine update_unevolved_species(state)
    ! STUB FOR INTEGRATOR
    type(burn_t)     :: state
    return
  end subroutine update_unevolved_species

  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    double precision :: Y(nspec)
    double precision :: raw_rates(4, nrates)
    double precision :: reactvec(num_rate_groups+2)
    integer :: i, j
    double precision :: dens, temp, rhoy

    Y(:) = state%xn(:) * aion_inv(:)
    dens = state%rho
    temp = state%T
    rhoy = dens*state%y_e
    
    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, temp, dens, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, temp, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(1:4)
    end do


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
         rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates

  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn
    use burn_type_module, only: net_itemp, net_ienuc

    implicit none

    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    type(plasma_state) :: pstate
    double precision :: Y(nspec)
    double precision :: ydot_nuc(nspec)
    double precision :: reactvec(num_rate_groups+2)
    integer :: i, j
    double precision :: dens, temp, rhoy, ye, enuc
    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)

    dens = state%rho
    temp = state%T

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(ydot_nuc, Y, rate_eval % screened_rates, dens)
    state%ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)
    
    ! weak Q-value modification dqweak (density and temperature dependent)
    
    ! weak particle energy generation rates from gamma heating and neutrino loss
    ! (does not include plasma neutrino losses)


    ! Get the neutrino losses
    call sneut5(temp, dens, state%abar, state%zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)
    state%ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state)
    
    ! write(*,*) '______________________________'
    ! do i = 1, nspec+2
    !    write(*,*) 'state%ydot(',i,'): ',state%ydot(i)
    ! end do
  end subroutine actual_rhs

  subroutine rhs_nuc(ydot_nuc, Y, screened_rates, dens)

    !$acc routine seq

    
    double precision, intent(out) :: ydot_nuc(nspec)
    double precision, intent(in)  :: Y(nspec)
    double precision, intent(in)  :: screened_rates(nrates)
    double precision, intent(in)  :: dens

    double precision :: scratch_0
    double precision :: scratch_1
    double precision :: scratch_2
    double precision :: scratch_3
    double precision :: scratch_4
    double precision :: scratch_5
    double precision :: scratch_6
    double precision :: scratch_7
    double precision :: scratch_8
    double precision :: scratch_9

    scratch_0 = screened_rates(k_n_p)*Y(jn)
    scratch_1 = Y(jc12)**2*dens
    scratch_2 = screened_rates(k_c12_c12n_mg23)*scratch_1
    scratch_3 = 0.5d0*scratch_2
    scratch_4 = screened_rates(k_c12_c12p_na23)*scratch_1
    scratch_5 = 0.5d0*scratch_4
    scratch_6 = screened_rates(k_c12_c12a_ne20)*scratch_1
    scratch_7 = 0.5d0*scratch_6
    scratch_8 = screened_rates(k_c12_ag_o16)*Y(jc12)*Y(jhe4)*dens
    scratch_9 = -scratch_8

    ydot_nuc(jn) = ( &
      -scratch_0 + scratch_3 &
       )

    ydot_nuc(jp) = ( &
      scratch_0 + scratch_5 &
       )

    ydot_nuc(jhe4) = ( &
      scratch_7 + scratch_9 &
       )

    ydot_nuc(jc12) = ( &
      -scratch_2 - scratch_4 - scratch_6 + scratch_9 &
       )

    ydot_nuc(jo16) = ( &
      scratch_8 &
       )

    ydot_nuc(jne20) = ( &
      scratch_7 &
       )

    ydot_nuc(jna23) = ( &
      scratch_5 &
       )

    ydot_nuc(jmg23) = ( &
      scratch_3 &
       )


  end subroutine rhs_nuc

  
  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: net_itemp, net_ienuc
    
    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    type(plasma_state) :: pstate
    double precision :: reactvec(num_rate_groups+2)
    double precision :: screened_rates_dt(nrates)
    double precision :: dfdy_nuc(nspec, nspec)
    double precision :: Y(nspec)
    double precision :: dens, temp, ye, rhoy, b1
    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer :: i, j

    dens = state%rho
    temp = state%T

    ! Set molar abundances
    Y(:) = state%xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)
    
    ! Species Jacobian elements with respect to other species
    call jac_nuc(dfdy_nuc, Y, rate_eval % screened_rates, dens)
    state%jac(1:nspec, 1:nspec) = dfdy_nuc

    ! Species Jacobian elements with respect to energy generation rate
    state%jac(1:nspec, net_ienuc) = 0.0d0

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
         rate_eval % unscreened_rates(i_dscor_dt, :) + &
         rate_eval % unscreened_rates(i_drate_dt, :) * &
         rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state%jac(1:nspec, net_itemp), Y, screened_rates_dt, dens)
    
    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses
    call sneut5(temp, dens, state%abar, state%zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    do j = 1, nspec
       b1 = ((aion(j) - state%abar) * state%abar * snuda + (zion(j) - state%zbar) * state%abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Energy generation rate Jacobian element with respect to energy generation rate
    state%jac(net_ienuc, net_ienuc) = 0.0d0

    ! Energy generation rate Jacobian element with respect to temperature
    call ener_gener_rate(state%jac(1:nspec, net_itemp), state%jac(net_ienuc, net_itemp))
    state%jac(net_ienuc, net_itemp) = state%jac(net_ienuc, net_itemp) - dsneutdt

    ! Add dqweak and epart contributions!!!

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac

  subroutine jac_nuc(dfdy_nuc, Y, screened_rates, dens)

    !$acc routine seq
    

    double precision, intent(out) :: dfdy_nuc(nspec, nspec)
    double precision, intent(in)  :: Y(nspec)
    double precision, intent(in)  :: screened_rates(nrates)
    double precision, intent(in)  :: dens

    double precision :: scratch_0
    double precision :: scratch_1
    double precision :: scratch_2
    double precision :: scratch_3
    double precision :: scratch_4
    double precision :: scratch_5
    double precision :: scratch_6
    double precision :: scratch_7
    double precision :: scratch_8
    double precision :: scratch_9

    scratch_0 = 1.0d0*Y(jc12)*dens
    scratch_1 = screened_rates(k_c12_c12n_mg23)*scratch_0
    scratch_2 = screened_rates(k_c12_c12p_na23)*scratch_0
    scratch_3 = screened_rates(k_c12_ag_o16)*dens
    scratch_4 = Y(jc12)*scratch_3
    scratch_5 = -scratch_4
    scratch_6 = Y(jhe4)*scratch_3
    scratch_7 = -scratch_6
    scratch_8 = screened_rates(k_c12_c12a_ne20)*scratch_0
    scratch_9 = 2.0d0*Y(jc12)*dens

    dfdy_nuc(jn,jn) = ( &
      -screened_rates(k_n_p) &
       )

    dfdy_nuc(jn,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn,jc12) = ( &
      scratch_1 &
       )

    dfdy_nuc(jn,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jn,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jn) = ( &
      screened_rates(k_n_p) &
       )

    dfdy_nuc(jp,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jc12) = ( &
      scratch_2 &
       )

    dfdy_nuc(jp,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jp,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jhe4) = ( &
      scratch_5 &
       )

    dfdy_nuc(jhe4,jc12) = ( &
      scratch_7 + scratch_8 &
       )

    dfdy_nuc(jhe4,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jhe4,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jhe4) = ( &
      scratch_5 &
       )

    dfdy_nuc(jc12,jc12) = ( &
      -screened_rates(k_c12_c12a_ne20)*scratch_9 - screened_rates(k_c12_c12n_mg23)*scratch_9 - &
      screened_rates(k_c12_c12p_na23)*scratch_9 + scratch_7 &
       )

    dfdy_nuc(jc12,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jc12,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jhe4) = ( &
      scratch_4 &
       )

    dfdy_nuc(jo16,jc12) = ( &
      scratch_6 &
       )

    dfdy_nuc(jo16,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jo16,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jc12) = ( &
      scratch_8 &
       )

    dfdy_nuc(jne20,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jne20,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jc12) = ( &
      scratch_2 &
       )

    dfdy_nuc(jna23,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jna23,jmg23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jn) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jp) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jhe4) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jc12) = ( &
      scratch_1 &
       )

    dfdy_nuc(jmg23,jo16) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jne20) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jna23) = ( &
      0.0d0 &
       )

    dfdy_nuc(jmg23,jmg23) = ( &
      0.0d0 &
       )

    
  end subroutine jac_nuc

end module actual_rhs_module
