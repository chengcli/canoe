// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// opacity
#include "freedman.hpp"

// coefficient from Richard S. Freedman 2014. APJS

Real c1 = 10.602, c2 = 2.882, c3 = 6.09e-15, c4 = 2.954, c5 = -2.526,
     c6 = 0.843, c7 = -5.490;
Real c8, c9, c10, c11, c12, c13 = 0.8321;

Real FreedmanMean2::GetAttenuation(Real wave1, Real wave2,
                                  AirParcel const& var) const {
  Real p = var.w[IPR];
  Real T = var.w[IDN];
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
	Real met = params_at("met");
	Real scale = params_at("scale");

  Real klowp = c1 * atan(logT - c2) -
               c3 / (logp + c4) * exp(pow(logT - c5, 2.0)) + c7;  // Eqn 4
	
	// xiz changes
	klowp += c6 * met;

  Real khigp = c8 + c9 * logT + c10 * pow(logT, 2.) +
               logp * (c11 + c12 * logT);  // Eqn 5
	
	// xiz changes
	khigp += c13*met*(0.5+1./M_PI*atan((logT-2.5)/0.2)); // Eqn 5

  Real result = pow(10.0, klowp) + pow(10.0, khigp);  // cm^2/g

  auto pthermo = Thermodynamics::GetInstance();
  Real dens = p / (pthermo->GetRd() * T);  // kg/m^3

  if (p > 5e1)
    return scale * 0.1 * dens * result;  // -> 1/m
  else
    return 0.;
}
