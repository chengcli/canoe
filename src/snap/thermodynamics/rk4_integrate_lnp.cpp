// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::rk4IntegrateLnp(AirParcel *air, Real dlnp, Method method,
                                     Real adlnTdlnP) const {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = air->w[IDN];
  Real pres = air->w[IPR];
  Real chi[4];
  Real latent[1 + NVAPOR];

  for (int rk = 0; rk < 4; ++rk) {
    EquilibrateTP(air);
    if (method != Method::ReversibleAdiabat) {
      for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
    }

    for (int i = 1; i <= NVAPOR; ++i) {
      // make a trial run to get latent heat
      air->w[i] += 1.E-6;
      auto rates = TryEquilibriumTP_VaporCloud(*air, i);
      latent[i] = GetLatentHeatMole(i, rates, air->w[IDN]) /
                  (Constants::Rgas * air->w[IDN]);
      air->w[i] -= 1.E-6;
    }

    // calculate tendency
    if (method == Method::ReversibleAdiabat ||
        method == Method::PseudoAdiabat) {
      chi[rk] = calDlnTDlnP(*air, latent);
    } else if (method == Method::DryAdiabat) {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = calDlnTDlnP(*air, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }
    chi[rk] += adlnTdlnP;

    // integrate over dlnp
    if (rk < 3) {
      air->w[IDN] = temp * exp(chi[rk] * dlnp * step[rk]);
    } else {
      air->w[IDN] =
          temp *
          exp(1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp);
    }
    air->w[IPR] = pres * exp(dlnp);
  }

  // recondensation
  EquilibrateTP(air);
  if (method != Method::ReversibleAdiabat) {
    for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
  }
}
