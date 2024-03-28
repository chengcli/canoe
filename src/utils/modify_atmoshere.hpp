#ifndef SRC_UTILS_MODIFY_ATMOSPHERE_HPP_
#define SRC_UTILS_MODIFY_ATMOSPHERE_HPP_

// C/C++
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>

// athena
#include <athena/mesh/mesh.hpp>

// helper functions, will be moved in the future
int find_pressure_level_lesser_pybind(Real pres, AthenaArray<Real> const &w, int k, int j, int is, int ie) ;

// modify atmoshere with adlnTdlnP
void modify_atmoshere_adlnTdlnP(MeshBlock *pmb, Real adlnTdlnP, Real pmin, Real pmax) ;

// modify atmoshere with adlnNH3dlnP
void modify_atmoshere_adlnNH3dlnP(MeshBlock *pmb, Real adlnNH3dlnP, Real pmin, Real pmax) ;

#endif //SRC_UTILS_MODIFY_ATMOSPHERE_HPP_

