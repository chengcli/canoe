// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::rk4IntegrateLnp(Variable *qfrac, Real dlnp, Method method,
                                     Real adlnTdlnP) const {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = qfrac->w[IDN];
  Real pres = qfrac->w[IPR];
  Real chi[4];
  Real latent[1 + NVAPOR];

  for (int rk = 0; rk < 4; ++rk) {
    EquilibrateTP(qfrac);
    if (method != Method::ReversibleAdiabat)
      for (int j = 0; j < NCLOUD; ++j) qfrac->c[j] = 0;

    // TODO(cli) : latent heat was diasabled now
    for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;

    // calculate tendency
    if (method == Method::ReversibleAdiabat ||
        method == Method::PseudoAdiabat) {
      chi[rk] = calDlnTDlnP(*qfrac, latent);
    } else if (method == Method::DryAdiabat) {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = calDlnTDlnP(*qfrac, latent);
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
  EquilibrateTP(qfrac);
  if (method != Method::ReversibleAdiabat)
    for (int j = 0; j < NCLOUD; ++j) qfrac->c[j] = 0;
}
