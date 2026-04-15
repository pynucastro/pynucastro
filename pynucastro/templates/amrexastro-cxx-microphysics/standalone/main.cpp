#include <iostream>

#include <extern_parameters.H>
#include <burn_type.H>
#include <eos.H>
#include <network.H>
#include <unit_test.H>
#include <actual_rhs.H>

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

    amrex::Array1D<amrex::Real, 1, neqs> ydot{};

    actual_rhs(burn_state, ydot);

    std::cout << std::setprecision(12);

    std::cout << "rho, T = " << burn_state.rho << " " << burn_state.T << std::endl;

    for (int n = 1; n <= NumSpec; ++n) {
        std::cout << "Ydot(" << short_spec_names_cxx[n-1] << ") = " << ydot(n) << std::endl;
    }

    std::cout << std::endl;



    amrex::Finalize();
}
