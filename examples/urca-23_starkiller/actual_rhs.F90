module actual_rhs_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use physical_constants, only: N_AVO
  use network
  use table_rates
  use burn_type_module

  implicit none

  type :: rate_eval_t
     real(rt) :: unscreened_rates(num_rate_groups, nrates)
     real(rt) :: screened_rates(nrates)
     real(rt) :: add_energy_rate(nrat_tabular)
  end type rate_eval_t
  
contains

  subroutine actual_rhs_init()
    ! STUB FOR MAESTRO'S TEST_REACT. ALL THE INIT IS DONE BY BURNER_INIT
    return
  end subroutine actual_rhs_init


  subroutine update_unevolved_species(state)
    ! STUB FOR INTEGRATOR
    type(burn_t)     :: state

    !$gpu
    
    return
  end subroutine update_unevolved_species


  subroutine zero_rate_eval(rate_eval)

    implicit none

    type(rate_eval_t), intent(inout) :: rate_eval

    !$gpu

    rate_eval % unscreened_rates(i_rate, :) = ZERO
    rate_eval % unscreened_rates(i_drate_dt, :) = ZERO
    rate_eval % unscreened_rates(i_scor, :) = ONE
    rate_eval % unscreened_rates(i_dscor_dt, :) = ZERO
    rate_eval % screened_rates = ZERO
    rate_eval % add_energy_rate = ZERO

  end subroutine zero_rate_eval


  subroutine evaluate_rates(state, rate_eval)
    !$acc routine seq

    use reaclib_rates, only: screen_reaclib, reaclib_evaluate
    use screening_module, only: screen5, plasma_state, fill_plasma_state

    implicit none
    
    type(burn_t)     :: state
    type(rate_eval_t), intent(out) :: rate_eval
    type(plasma_state) :: pstate
    real(rt) :: Y(nspec)
    real(rt) :: reactvec(num_rate_groups)
    integer :: i, j
    real(rt) :: rhoy, scor, dscor_dt, dscor_dd

    !$gpu

    Y(:) = state % xn(:) * aion_inv(:)
    rhoy = state % rho * state % y_e

    ! Zero out the rates
    call zero_rate_eval(rate_eval)

    ! Calculate Reaclib rates
    call fill_plasma_state(pstate, state % T, state % rho, Y)
    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, state % T, i, reactvec)
       rate_eval % unscreened_rates(:,i) = reactvec(:)
    end do

    ! Evaluate screening factors
    if (screen_reaclib) then

      call screen5(pstate, 1, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,1) = scor
      rate_eval % unscreened_rates(i_dscor_dt,1) = dscor_dt
      rate_eval % unscreened_rates(i_scor,2) = scor
      rate_eval % unscreened_rates(i_dscor_dt,2) = dscor_dt
      rate_eval % unscreened_rates(i_scor,3) = scor
      rate_eval % unscreened_rates(i_dscor_dt,3) = dscor_dt


      call screen5(pstate, 2, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,4) = scor
      rate_eval % unscreened_rates(i_dscor_dt,4) = dscor_dt

    end if

    ! Calculate tabular rates
    call tabular_evaluate(rate_table_j_na23_ne23, rhoy_table_j_na23_ne23, temp_table_j_na23_ne23, &
                          num_rhoy_j_na23_ne23, num_temp_j_na23_ne23, num_vars_j_na23_ne23, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(i_rate:i_scor,6) = reactvec(i_rate:i_scor)
    rate_eval % add_energy_rate(1)  = reactvec(i_eneut)

    call tabular_evaluate(rate_table_j_ne23_na23, rhoy_table_j_ne23_na23, temp_table_j_ne23_na23, &
                          num_rhoy_j_ne23_na23, num_temp_j_ne23_na23, num_vars_j_ne23_na23, &
                          rhoy, state % T, reactvec)
    rate_eval % unscreened_rates(i_rate:i_scor,7) = reactvec(i_rate:i_scor)
    rate_eval % add_energy_rate(2)  = reactvec(i_eneut)


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state)
    
    !$acc routine seq

    use extern_probin_module, only: do_constant_volume_burn, disable_thermal_neutrinos
    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: Y(nspec), ydot_nuc(nspec)
    real(rt) :: reactvec(num_rate_groups)
    integer :: i, j
    real(rt) :: rhoy, ye, enuc
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(state, ydot_nuc, Y, rate_eval % screened_rates)
    state % ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)

    ! include reaction neutrino losses (non-thermal)
    enuc = enuc + N_AVO * Y(jna23) * rate_eval % add_energy_rate(j_na23_ne23)
    enuc = enuc + N_AVO * Y(jne23) * rate_eval % add_energy_rate(j_ne23_na23)

    ! Get the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    else
       sneut = ZERO
    end if

    ! Append the energy equation (this is erg/g/s)
    state % ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state)

  end subroutine actual_rhs


  subroutine rhs_nuc(state, ydot_nuc, Y, screened_rates)

    !$acc routine seq

    implicit none

    type (burn_t), intent(in) :: state
    real(rt), intent(out) :: ydot_nuc(nspec)
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)

    !$gpu



    ydot_nuc(jn) = ( &
      0.5d0*screened_rates(k_c12_c12__n_mg23)*Y(jc12)**2*state % rho - &
      screened_rates(k_n__p__weak__wc12)*Y(jn) &
       )

    ydot_nuc(jp) = ( &
      0.5d0*screened_rates(k_c12_c12__p_na23)*Y(jc12)**2*state % rho + &
      screened_rates(k_n__p__weak__wc12)*Y(jn) &
       )

    ydot_nuc(jhe4) = ( &
      0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho - &
      screened_rates(k_c12_c12__n_mg23)*Y(jc12)**2*state % rho - &
      screened_rates(k_c12_c12__p_na23)*Y(jc12)**2*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )

    ydot_nuc(jne20) = ( &
      0.5d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)**2*state % rho &
       )

    ydot_nuc(jne23) = ( &
      screened_rates(k_na23__ne23)*Y(jna23) - screened_rates(k_ne23__na23)*Y(jne23) &
       )

    ydot_nuc(jna23) = ( &
      0.5d0*screened_rates(k_c12_c12__p_na23)*Y(jc12)**2*state % rho - &
      screened_rates(k_na23__ne23)*Y(jna23) + screened_rates(k_ne23__na23)*Y(jne23) &
       )

    ydot_nuc(jmg23) = ( &
      0.5d0*screened_rates(k_c12_c12__n_mg23)*Y(jc12)**2*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state)

    !$acc routine seq

    use extern_probin_module, only: disable_thermal_neutrinos
    use burn_type_module, only: net_itemp, net_ienuc
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t) :: state
    type(rate_eval_t) :: rate_eval
    real(rt) :: screened_rates_dt(nrates)
    real(rt) :: Y(nspec), yderivs(nspec)
    real(rt) :: ye, rhoy, b1, scratch
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer  :: j, k

    !$gpu

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)

    ! Zero out the Jacobian
    call set_jac_zero(state)

    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, Y, rate_eval % screened_rates)

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
                        rate_eval % unscreened_rates(i_dscor_dt, :) + &
                        rate_eval % unscreened_rates(i_drate_dt, :) * &
                        rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, yderivs, Y, screened_rates_dt)

    do k = 1, nspec
       call set_jac_entry(state, k, net_itemp, yderivs(k))
    enddo

    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(state, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

       do j = 1, nspec
          b1 = ((aion(j) - state % abar) * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
          call get_jac_entry(state, net_ienuc, j, scratch)
          scratch = scratch - b1
          call set_jac_entry(state, net_ienuc, j, scratch)
       enddo
    endif

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(state, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    if (.not. disable_thermal_neutrinos) then
       scratch = scratch - dsneutdt
    endif
    call set_jac_entry(state, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac


  subroutine jac_nuc(state, Y, screened_rates)

    !$acc routine seq

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(inout) :: state
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch


    !$gpu


    scratch = (&
      -screened_rates(k_n__p__weak__wc12) &
       )
    call set_jac_entry(state, jn, jn, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__n_mg23)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jn, jc12, scratch)

    scratch = (&
      screened_rates(k_n__p__weak__wc12) &
       )
    call set_jac_entry(state, jp, jn, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__p_na23)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jp, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jhe4, jhe4, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jc12, jhe4, scratch)

    scratch = (&
      -2.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*state % rho - 2.0d0* &
      screened_rates(k_c12_c12__n_mg23)*Y(jc12)*state % rho - 2.0d0* &
      screened_rates(k_c12_c12__p_na23)*Y(jc12)*state % rho - &
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jo16, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(state, jo16, jc12, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__he4_ne20)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jne20, jc12, scratch)

    scratch = (&
      -screened_rates(k_ne23__na23) &
       )
    call set_jac_entry(state, jne23, jne23, scratch)

    scratch = (&
      screened_rates(k_na23__ne23) &
       )
    call set_jac_entry(state, jne23, jna23, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__p_na23)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jna23, jc12, scratch)

    scratch = (&
      screened_rates(k_ne23__na23) &
       )
    call set_jac_entry(state, jna23, jne23, scratch)

    scratch = (&
      -screened_rates(k_na23__ne23) &
       )
    call set_jac_entry(state, jna23, jna23, scratch)

    scratch = (&
      1.0d0*screened_rates(k_c12_c12__n_mg23)*Y(jc12)*state % rho &
       )
    call set_jac_entry(state, jmg23, jc12, scratch)


  end subroutine jac_nuc

end module actual_rhs_module
