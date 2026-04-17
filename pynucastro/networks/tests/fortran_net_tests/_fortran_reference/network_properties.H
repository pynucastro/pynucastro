#ifndef NETWORK_PROPERTIES_H
#define NETWORK_PROPERTIES_H
#include <vector>
#include <string>
#include <amrex_bridge.H>

constexpr int NumSpec = 9;

// Note: these are 0-based
constexpr Real aion[NumSpec] = {
       1.0, // 0 : n
       1.0, // 1 : p
       4.0, // 2 : He4
      12.0, // 3 : C12
      16.0, // 4 : O16
      20.0, // 5 : Ne20
      23.0, // 6 : Na23
      23.0, // 7 : Mg23
      24.0, // 8 : Mg24
 };

constexpr Real aion_inv[NumSpec] = {
    1.0/   1.0, // 0 : n
    1.0/   1.0, // 1 : p
    1.0/   4.0, // 2 : He4
    1.0/  12.0, // 3 : C12
    1.0/  16.0, // 4 : O16
    1.0/  20.0, // 5 : Ne20
    1.0/  23.0, // 6 : Na23
    1.0/  23.0, // 7 : Mg23
    1.0/  24.0, // 8 : Mg24
 };

constexpr Real zion[NumSpec] = {
       0.0, // 0 : n
       1.0, // 1 : p
       2.0, // 2 : He4
       6.0, // 3 : C12
       8.0, // 4 : O16
      10.0, // 5 : Ne20
      11.0, // 6 : Na23
      12.0, // 7 : Mg23
      12.0, // 8 : Mg24
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
    "Mg24", // 8
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
    Mg24,
  };
}

#endif
