#ifndef SRC_UTILS_CONSTRUCT_ATMOSPHERE_DRY_HPP_
#define SRC_UTILS_CONSTRUCT_ATMOSPHERE_DRY_HPP_

// C/C++
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>


// set up an adiabatic atmosphere
void construct_atmosphere(MeshBlock *pmb, ParameterInput *pin, Real xNH3, Real T0);

#endif  // SRC_UTILS_CONSTRUCT_ATMOSPHERE_HPP_