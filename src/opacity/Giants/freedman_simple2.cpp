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
#include "freedman.hpp"

// xiz semigrey
Real FreedmanSimple2::GetAttenuation(Real wave1, Real wave2,
                                     AirParcel const &var) const {
  auto pthermo = Thermodynamics::GetInstance();
  Real mu = Constants::Rgas / pthermo->GetRd();
  Real result;
  Real p = var.w[IPR];
  Real T = var.w[IDN];
  Real scale = GetPar<Real>("scale");

  /*xiz semigrey
  //Tan and Komacek 2019 simple fit m^2/kg
  //    k_VIS= 10.0_dp**(0.0478_dp*Pl10**2 - 0.1366_dp*Pl10 - 3.2095_dp)
  //    k_IR = 10.0_dp**(0.0498_dp*Pl10**2 - 0.1329_dp*Pl10 - 2.9457_dp)
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
  // result = 1.e-6 * pow(p, 0.5);  // visible opacity scale in disort
  
  
  //Shami intercomparison 2024 from Guillot
    result = 1.e-3;

  Real dens = p * mu / (Constants::Rgas * T);  // kg/m^3
                                               //  if (p > 5e1)
  return scale * dens * result;                // -> 1/m
  //  else
  //    return 0.;
}
