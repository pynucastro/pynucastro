#include <amrex_bridge.H>
#include <network_properties.H>
#include <actual_rhs.H>

#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {

    // we can do
    //   ./main rho T
    // to pass in the density and temperature.
    // composition will always be uniform.

    // defaults
    double T{1.e9};
    double rho{2.e8};

    if (argc == 3) {
        rho = std::atof(argv[1]);
        T = std::atof(argv[2]);
    } else if (argc != 1) {
        std::cout << "Error: invalid inputs" << std::endl;
        std::exit(1);
    }

    std::cout << "(rho, T) = " << rho << " " << T << std::endl;

    actual_network_init();

    burn_t state;
    state.T = T;
    state.rho = rho;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = 1.0 / static_cast<Real>(NumSpec);
    }

    Array1D<Real, 1, NumSpec> ydot;
    actual_rhs(state, ydot);

    std::cout << std::setprecision(12);

    for (int n = 1; n <= NumSpec; ++n) {
        std::cout << "Ydot(" << spec_names[n-1] << ") = " << ydot(n) << std::endl;
    }

    std::cout << std::endl;

    MathArray2D<1, NumSpec, 1, NumSpec> jac;
    actual_jac(state, jac);

    for (int jcol = 1; jcol <= NumSpec; ++jcol) {
        for (int irow = 1; irow <= NumSpec; ++irow) {
            std::cout << "jac(" << irow << "," << jcol << ") = "
                      << jac(irow,jcol) << std::endl;
        }
    }

    std::cout << std::endl;

    Real enuc;
    ener_gener_rate(ydot, enuc);
    std::cout << "instantaneous energy generation rate (erg/g/s) = " << enuc << std::endl;

}
