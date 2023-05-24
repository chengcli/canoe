// C/C++ headers
#include "freedman_mean.hpp"

#include <algorithm>
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>
#include <cassert>  // assert
#include <cmath>
#include <cstring>
#include <iostream>
#include <snap/cell_variables.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>
#include <stdexcept>

// coefficient from Richard S. Freedman 2014. APJS

Real c1 = 10.602, c2 = 2.882, c3 = 6.09e-15, c4 = 2.954, c5 = -2.526,
     c6 = 0.843, c7 = -5.490;
Real c8, c9, c10, c11, c12;

Real FreedmanMean::getAttenuation(Real wave1, Real wave2,
                                  CellVariables const& var) const {
  Real p = var.q[IPR];
  Real T = var.q[IDN];
  if (T < 800.) {
    c8 = -14.051;
    c9 = 3.055;
    c10 = 0.024;
    c11 = 1.877;
    c12 = -0.445;
  } else {
    c8 = 82.241;
    c9 = -55.456;
    c10 = 8.754;
    c11 = 0.7048;
    c12 = -0.0414;
  }
  Real logp = log10(p * 10.);  // Pa to dyn/cm2
  Real logT = log10(T);

  Real klowp = c1 * atan(logT - c2) -
               c3 / (logp + c4) * exp(pow(logT - c5, 2.0)) + c7;  // Eqn 4

  Real khigp = c8 + c9 * logT + c10 * pow(logT, 2.) +
               logp * (c11 + c12 * logT);  // Eqn 5

  Real result = pow(10.0, klowp) + pow(10.0, khigp);  // cm^2/g

  Real dens = p / (pthermo_->GetRd() * T);  // kg/m^3

  if (p > 5e1)
    return 0.1 * dens * result;  // -> 1/m
  else
    return 0.;
}
