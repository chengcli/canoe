// C/C++
#include <algorithm>
#include <iomanip>
#include <iostream>

// thermodynamics
#include "atm_thermodynamics.hpp"

void rk4_integrate_z(AirParcel* air, Real dz, std::string method, Real grav,
                     Real adTdz) {
  auto pthermo = Thermodynamics::GetInstance();

  auto const& mu_ratio = pthermo->GetMuRatio();
  auto const& cp_ratio_mole = pthermo->GetCpRatioMole();

  Real step[] = {0.5, 0.5, 1.};
  Real temp = air->w[IDN];
  Real pres = air->w[IPR];
  Real dTdz[4], chi[4];
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

    Real q_gas = 1., q_eps = 1.;
    for (int i = 1; i <= NVAPOR; ++i) {
      q_eps += air->w[i] * (mu_ratio[i] - 1.);
    }

    for (int j = 0; j < NCLOUD; ++j) {
      q_eps += air->c[j] * (mu_ratio[1 + NVAPOR + j] - 1.);
      q_gas += -air->c[j];
    }

    Real g_ov_Rd = grav / pthermo->GetRd();
    Real R_ov_Rd = q_gas / q_eps;

    if (method == "reversible" || method == "pseudo") {
      chi[rk] = cal_dlnT_dlnP(*air, cp_ratio_mole, latent);
    } else if (method == "dry") {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = cal_dlnT_dlnP(*air, cp_ratio_mole, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }

    dTdz[rk] = -chi[rk] * g_ov_Rd / R_ov_Rd + adTdz;
    chi[rk] = -R_ov_Rd / g_ov_Rd * dTdz[rk];

    // integrate over dz
    Real chi_avg;
    if (rk < 3) {
      air->w[IDN] = temp + dTdz[rk] * dz * step[rk];
      chi_avg = chi[rk];
    } else {
      air->w[IDN] =
          temp +
          1. / 6. * (dTdz[0] + 2. * dTdz[1] + 2. * dTdz[2] + dTdz[3]) * dz;
      chi_avg = 1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]);
    }

    if (!(air->w[IDN] > 0.)) air->w[IDN] = temp;

    if (fabs(air->w[IDN] - temp) > 0.01) {
      air->w[IPR] = pres * pow(air->w[IDN] / temp, 1. / chi_avg);
    } else {  // isothermal limit
      air->w[IPR] =
          pres * exp(-2. * g_ov_Rd * dz / (R_ov_Rd * (air->w[IDN] + temp)));
    }
  }

  // recondensation
  pthermo->EquilibrateTP(air);
  if (method != "reversible") {
    for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
  }
}
