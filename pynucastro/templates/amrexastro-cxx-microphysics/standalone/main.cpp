#include <iostream>

#include <extern_parameters.H>
#include <burn_type.H>
#include <eos.H>
#include <network.H>
#include <unit_test.H>
#include <actual_rhs.H>
#include <ArrayUtilities.H>

int main(int argc, char *argv[]) {

    amrex::Initialize(argc, argv);

    std::cout << "starting the single zone burn..." << std::endl;

    init_unit_test();

    // C++ EOS initialization (must be done after Fortran eos_init and
    // init_extern_parameters)
    eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);

    // C++ Network, RHS, screening, rates initialization
    network_init();

    // setup the state for the burn

    burn_t burn_state;

    burn_state.rho = testing_rp::density;
    burn_state.T = testing_rp::temperature;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = 1.0 / NumSpec;
    }

    // call the EOS -- this will set the composition variables
    // (including Ye) and the specific heat (c_v)
    eos(eos_input_rt, burn_state);

    // get the Ydots

    amrex::Array1D<amrex::Real, 1, neqs> ydot{};

    actual_rhs(burn_state, ydot);

    std::cout << std::setprecision(14);

    std::cout << std::endl;
    std::cout << "rho, T = " << burn_state.rho << " " << burn_state.T << std::endl;
    std::cout << std::endl;

    for (int n = 1; n <= NumSpec; ++n) {
        std::cout << "Ydot(" << short_spec_names_cxx[n-1] << ") = "
                  << std::setw(5) << ydot(n) << std::endl;
    }

    std::cout << std::endl;

    // get the Jacobian
    // note that this is in terms of Y and e

    ArrayUtil::MathArray2D<amrex::Real, 1, NumSpec+1, 1, NumSpec+1> jac;
    jac.zero();

    actual_jac(burn_state, jac);

    std::cout << "Jacobian values" << std::endl;

    for (int irow = 1; irow <= NumSpec+1; ++irow) {
        auto irow_name = irow <= NumSpec ? short_spec_names_cxx[irow-1] : "e";

        for (int jcol = 1; jcol <= NumSpec+1; ++jcol) {
            auto jcol_name = jcol <= NumSpec ? short_spec_names_cxx[jcol-1] : "e";
            std::cout << "jac("
                      << std::setw(5) << irow_name << ","
                      << std::setw(5) << jcol_name << ") = "
                      << jac(irow, jcol) << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    // get the rates -- this is just the N_A<σv>

    // create molar fractions
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = burn_state.xn[n-1] * aion_inv[n-1];
    }

    // compute and output energy generation rates

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives{1};
    evaluate_rates<do_T_derivatives>(burn_state, Y, rate_eval);

    amrex::Real enuc{};
    ener_gener_rate(ydot, enuc);

    std::cout << "Instantaneous energy generation rate" << std::endl;
    std::cout << "ε_nuc = " << enuc << std::endl;
    std::cout << "ε_{ν,weak} = " << rate_eval.enuc_weak << std::endl;
    std::cout << std::endl;

    // output reaction rates

    std::cout << "rates" << std::endl;
    for (int n = 1; n <= Rates::NumRates; ++n) {
        std::cout << "rate(" << std::setw(40) << Rates::rate_names[n] << ") = "
                  << rate_eval.screened_rates(n) << std::endl;
    }

    std::cout << std::endl;

    amrex::Finalize();
}
