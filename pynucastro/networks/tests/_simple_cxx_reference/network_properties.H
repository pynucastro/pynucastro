#ifndef NETWORK_PROPERTIES_H
#define NETWORK_PROPERTIES_H
#include <vector>
#include <string>
#include <amrex_bridge.H>

constexpr int NumSpec = 8;

constexpr Real aion[NumSpec] = {
       1.0, // 0
       1.0, // 1
       4.0, // 2
      12.0, // 3
      16.0, // 4
      20.0, // 5
      23.0, // 6
      23.0, // 7
 };

constexpr Real aion_inv[NumSpec] = {
    1.0/   1.0, // 0
    1.0/   1.0, // 1
    1.0/   4.0, // 2
    1.0/  12.0, // 3
    1.0/  16.0, // 4
    1.0/  20.0, // 5
    1.0/  23.0, // 6
    1.0/  23.0, // 7
 };

constexpr Real zion[NumSpec] = {
       0.0, // 0
       1.0, // 1
       2.0, // 2
       6.0, // 3
       8.0, // 4
      10.0, // 5
      11.0, // 6
      12.0, // 7
 };

static const std::vector<std::string> spec_names = {
    "N", // 0
    "H1", // 1
    "He4", // 2
    "C12", // 3
    "O16", // 4
    "Ne20", // 5
    "Na23", // 6
    "Mg23", // 7
 };

namespace Species {
  enum NetworkSpecies {
    N=1,
    H1,
    He4,
    C12,
    O16,
    Ne20,
    Na23,
    Mg23,
  };
}

#endif
