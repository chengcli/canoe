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
#include <impl.hpp>
#include <variable.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include "freedman_simple.hpp"

FreedmanSimple::FreedmanSimple(MeshBlock *pmb, ParameterInput *pin,
                               std::string bname)
    : Absorber("FreedmanSimple") {
  char str[80];
  snprintf(str, sizeof(str), "%s.%s.scale", bname.c_str(), name_.c_str());
  scale_ = pin->GetOrAddReal("radiation", str, 1.);
}

// xiz semigrey
Real FreedmanSimple::GetAttenuation(Real wave1, Real wave2,
                                    Variable const &var) const {
  static const Real Rgas = 8.314462;
  const Absorber *pabs = this;
  Real mu = Rgas / pmy_block_->pimpl->pthermo->GetRd();
  Real result;
  Real p = var.q[IPR];
  Real T = var.q[IDN];

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
  result = 2.28e-6 * pow(p, 0.53);  // visible opacity scale in disort
  //

  Real dens = p * mu / (Rgas * T);  // kg/m^3
                                    //  if (p > 5e1)
  return scale_ * dens * result;    // -> 1/m
  //  else
  //    return 0.;
}
