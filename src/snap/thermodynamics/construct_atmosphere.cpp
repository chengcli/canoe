// C/C++
#include <cstdlib>
#include <iostream>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::ConstructAtmosphere(Variable *qfrac, Real dzORdlnp,
                                         Method method, Real userp) const {
  Real gammad = pmy_block_->peos->GetGamma();
  Real grav = -pmy_block_->phydro->hsrc.GetG1();
  auto pimpl = pmy_block_->pimpl;

  int isat[1 + NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n) isat[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  Real xg = 1.;  // total moles in gas phase
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    Real rate = TryEquilibriumTP(*qfrac, iv);
    if (rate > 0.) isat[iv] = 1;
    ExecuteEquilibriumTP(qfrac, iv);
    xg -= rate;
  }

  // RK4 integration
  Real adlnTdlnP, adTdz;
#ifdef HYDROSTATIC
  adlnTdlnP = userp;
  rk4IntegrateLnp(qfrac, isat, gammad, dzORdlnp, method, adlnTdlnP);
#else
  adTdz = userp;
  rk4IntegrateZ(qfrac, isat, gammad, dzORdlnp, method, adTdz);
#endif
}
