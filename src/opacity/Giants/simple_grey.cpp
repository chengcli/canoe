// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// opacity
#include "grey_gas.hpp"

// xiz semigrey
Real SimpleGrey::GetAttenuation(Real wave1, Real wave2,
                                AirParcel const &var) const {
  auto pthermo = Thermodynamics::GetInstance();
  Real mu = Constants::Rgas / pthermo->GetRd();
  Real p = var.w[IPR];
  Real T = var.w[IDN];
  Real dens = p * mu / (Constants::Rgas * T);  // kg/m^3

  Real kappa_a = GetPar<Real>("kappa_a");
  Real kappa_b = GetPar<Real>("kappa_b");
  Real kappa_cut = GetPar<Real>("kappa_cut");

  // kappa = kappa_a * pow(p, kappa_b)
  Real kappa = kappa_a * pow(p, kappa_b);
  if (kappa < kappa_cut) {
    kappa = kappa_cut;
  }
  // Tan and Komacek 2019 simple fit m^2/kg
  //     k_VIS= 10.0_dp**(0.0478_dp*Pl10**2 - 0.1366_dp*Pl10 - 3.2095_dp)
  //     k_IR = 10.0_dp**(0.0498_dp*Pl10**2 - 0.1329_dp*Pl10 - 2.9457_dp)
  /*
    Real logp = log10(p); // Pa
  if (wave < 40000.) //for semigrey
    result = pow(10.0, (0.0498*pow(logp,2.) - 0.1329*logp - 2.9457));
  else
    result = pow(10.0, (0.0478*pow(logp,2.) - 0.1366*logp - 3.2095));
  */

  // Komacek et al. 2017
  // if (wave < 40000.) //for semigrey
  //   result = 2.28e-6*pow(p, 0.53);
  // else
  // result = 2.28e-6 * pow(p, 0.53);  // visible opacity scale in disort
  //

  // Shami intercomparison 2024 from Guillot
  // result = 1.e-3;

  return dens * kappa;  // -> 1/m
}
