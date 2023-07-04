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
  int isat[1 + NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n) isat[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  Real xg = 1.;  // total moles in gas phase
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    auto rates = TryEquilibriumTP(*qfrac, iv);

    // saturation indicator
    if (rates[0] < 0.) isat[iv] = 1;

    // vapor condensation rate
    qfrac->w[iv] += rates[0];

    // cloud concentration rates
    for (int n = 1; n < rates.size(); ++n)
      qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
    xg -= rate;
  }

  // RK4 integration
#ifdef HYDROSTATIC
  rk4IntegrateLnp(qfrac, isat, dzORdlnp, method, userp);
#else
  rk4IntegrateZ(qfrac, isat, dzORdlnp, method, userp);
#endif
}
