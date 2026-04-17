#ifndef BURN_TYPE_H
#define BURN_TYPE_H

#include <amrex_bridge.H>
#include <network_properties.H>

struct burn_t {

  Real rho;
  Real T;
  Real xn[NumSpec];
  Real y_e;

};

inline
void compute_ye(burn_t& state) {
    state.y_e = 0.0;
    for (int n = 0; n < NumSpec; ++n) {
        state.y_e += state.xn[n] * zion[n] * aion_inv[n];
    }
}

#endif
