#include <amrex_bridge.H>
#include <iostream>

#include <burn_type.H>
#include <actual_network.H>
#include <actual_rhs.H>

extern "C" {

    void
    network_init() {
        actual_network_init();
    }

    void
    rhs_f(Real rho, Real T, Real* X, Real* dYdt) {

        burn_t state;
        state.rho = rho;
        state.T = T;
        for (int n = 0; n < NumSpec; ++n) {
            state.xn[n] = X[n];
        }

        Array1D<Real, 1, NumSpec> ydot;

        actual_rhs(state, ydot);

        // copy back to 0-index dydt
        for (int n = 0; n < NumSpec; ++n) {
            dYdt[n] = ydot(n+1);
        }

    }


    void
    jac_f(Real rho, Real T, Real* X, Real* J) {

        MathArray2D<1, NumSpec, 1, NumSpec> jac;

        burn_t state;
        state.rho = rho;
        state.T = T;
        for (int n = 0; n < NumSpec; ++n) {
            state.xn[n] = X[n];
        }

        actual_jac(state, jac);

        // copy back to 0-indexed J in column-major order to be compatible
        // with Fortran

        for (int jcol = 0; jcol < NumSpec; ++jcol) {
            for (int irow = 0; irow < NumSpec; ++ irow) {
                J[jcol*NumSpec + irow] = jac(irow+1, jcol+1);
            }
        }

    }


    void
    ener_gener_f(Real* dYdt, Real& enuc) {

        Array1D<Real, 1, NumSpec> ydot;
        for (int n = 0; n < NumSpec; ++n) {
            ydot(n+1) = dYdt[n];
        }

        ener_gener_rate(ydot, enuc);

    }

}
