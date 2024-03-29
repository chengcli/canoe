// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "atm_thermodynamics.hpp"

void rk4_integrate_lnp(AirParcel* air, Real dlnp, std::string method,
                       Real adlnTdlnP) {
  auto pthermo = Thermodynamics::GetInstance();
  auto const& cp_ratio_mole = pthermo->GetCpRatioMole();

  Real step[] = {0.5, 0.5, 1.};
  Real temp = air->w[IDN];
  Real pres = air->w[IPR];
  Real chi[4];
  Real latent[1 + NVAPOR];

  for (int rk = 0; rk < 4; ++rk) {
    pthermo->EquilibrateTP(air);
    if (method != "reversible") {
      for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
    }

    for (int i = 1; i <= NVAPOR; ++i) {
      // make a trial run to get latent heat
      air->w[i] += 1.E-6;
      auto rates = pthermo->TryEquilibriumTP_VaporCloud(*air, i);
      latent[i] = pthermo->GetLatentHeatMole(i, rates, air->w[IDN]) /
                  (Constants::Rgas * air->w[IDN]);
      air->w[i] -= 1.E-6;
    }

    // calculate tendency
    if (method == "reversible" || method == "pseudo") {
      chi[rk] = cal_dlnT_dlnP(*air, cp_ratio_mole, latent);
    } else if (method == "dry") {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = cal_dlnT_dlnP(*air, cp_ratio_mole, latent);
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
  pthermo->EquilibrateTP(air);
  if (method != "reversible") {
    for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
  }
}
