// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"
#include "thermodynamics_helper.hpp"

void rk4IntegrateLnp(Variable *qfrac, int isat[] Real gamma, Real dlnp,
                     int method, Real adlnTdlnP) {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = qfrac.w[IDN];
  Real pres = qfrac.w[IPR];
  Real chi[4];

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor and cloud
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      setTotalEquivalentVapor(qfrac);
      isat[iv] = 0;
    }

    for (int iv = 1; iv <= NVAPOR; ++iv) {
      Real rate = TryEquilibriumTP(qfrac, iv);
      if (rate > 0.) isat[iv] = 1;
      if (method == ReversibleAdiabat) {
        ExecuteEquilibriumTP(qfrac);
      } else {
        qfrac.w[iv] -= rate;
      }
    }

    // calculate gamma
    gammad = updateGammad(qfrac);

    // calculate tendency
    if (method == ReversibleAdiabat || method == PseudoAdiabat)
      chi[rk] = CaldlnTdlnP(qfrac, isat, gammad);
    else if (method == DryAdiabat) {
      for (int iv = 1; iv <= NVAPOR; ++iv) isat[iv] = 0;
      chi[rk] = CaldlnTdlnP(qfrac, isat, gammad);
    } else {  // isothermal
      chi[rk] = 0.;
    }
    chi[rk] += adlnTdlnP;

    // integrate over dlnp
    if (rk < 3) {
      qfrac.w[IDN] = temp * exp(chi[rk] * dlnp * step[rk]);
    } else {
      qfrac.w[IDN] =
          temp *
          exp(1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp);
    }
    qfrac.w[IPR] = pres * exp(dlnp);
  }

  // recondensation
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    Real rate = TryEquilibriumTP(qfrac, iv);
    if (rate > 0.) isat[iv] = 1;
    if (method == ReversibleAdiabat) {
      ExecuteEquilibriumTP(qfrac);
    } else {
      qfrac.w[iv] -= rate;
    }
  }
}
