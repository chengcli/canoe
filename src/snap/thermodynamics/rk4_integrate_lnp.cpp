// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"
#include "thermodynamics_helper.hpp"

void Thermodynamics::rk4IntegrateLnp(Variable *qfrac, int isat[], Real dlnp,
                                     Method method, Real adlnTdlnP) {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = qfrac->w[IDN];
  Real pres = qfrac->w[IPR];
  Real chi[4];

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor and cloud
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      setTotalEquivalentVapor(qfrac, iv);
      isat[iv] = 0;
    }

    for (int iv = 1; iv <= NVAPOR; ++iv) {
      auto rates = TryEquilibriumTP(*qfrac, iv);

      // saturation indicator
      if (rates[0] > 0.) isat[iv] = 1;

      // vapor condensation rate
      qfrac->w[iv] -= rates[0];

      // cloud concentration rates
      if (method == ReversibleAdiabat) {
        for (int n = 1; n < rates.size(); ++n)
          qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
      }
    }

    // calculate tendency
    if (method == ReversibleAdiabat || method == PseudoAdiabat)
      chi[rk] = CaldlnTdlnP(*qfrac, isat);
    else if (method == DryAdiabat) {
      for (int iv = 1; iv <= NVAPOR; ++iv) isat[iv] = 0;
      chi[rk] = CaldlnTdlnP(*qfrac, isat);
    } else {  // isothermal
      chi[rk] = 0.;
    }
    chi[rk] += adlnTdlnP;

    // integrate over dlnp
    if (rk < 3) {
      qfrac->w[IDN] = temp * exp(chi[rk] * dlnp * step[rk]);
    } else {
      qfrac->w[IDN] =
          temp *
          exp(1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp);
    }
    qfrac->w[IPR] = pres * exp(dlnp);
  }

  // recondensation
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    auto rates = TryEquilibriumTP(*qfrac, iv);

    // saturation indicator
    if (rates[0] > 0.) isat[iv] = 1;

    // vapor condensation rate
    qfrac->w[iv] -= rates[0];

    // cloud concentration rates
    if (method == ReversibleAdiabat) {
      for (int n = 1; n < rates.size(); ++n)
        qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
    }
  }
}
