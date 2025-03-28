#ifndef BURN_TYPE_H
#define BURN_TYPE_H

#include <amrex_bridge.H>
#include <network_properties.H>

struct burn_t {

  Real rho;
  Real T;
  Real xn[NumSpec];

};

#endif
