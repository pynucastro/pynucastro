#include <amrex_bridge.H>
#include <iostream>
#include <cassert>

#include <burn_type.H>
#include <actual_rhs.H>

extern "C" {

    void
    rhs_f(Real rho, Real T, Real* X, Real* dYdt, int size) {

        assert(size == NumSpec);

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
    jac_f(Real rho, Real T, Real* X, Real* J, int size) {

        assert(size == NumSpec);

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
}
