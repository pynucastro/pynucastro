#include <amrex_bridge.H>
#include <network_properties.H>
#include <actual_rhs.H>

#include <iostream>
#include <iomanip>

int main() {

    actual_network_init();

    burn_t state;
    state.T = 1.e9;
    state.rho = 2.e8;
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
