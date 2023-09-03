// C/C++
#include <algorithm>
#include <iomanip>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::rk4IntegrateZ(AirParcel *air, Real dz, Method method,
                                   Real grav, Real adTdz) const {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = air->w[IDN];
  Real pres = air->w[IPR];
  Real dTdz[4], chi[4];
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

    Real q_gas = 1., q_eps = 1.;
    for (int i = 1; i <= NVAPOR; ++i) {
      q_eps += air->w[i] * (mu_ratio_[i] - 1.);
    }

    for (int j = 0; j < NCLOUD; ++j) {
      q_eps += air->c[j] * (mu_ratio_[1 + NVAPOR + j] - 1.);
      q_gas += -air->c[j];
    }

    Real g_ov_Rd = grav / Rd_;
    Real R_ov_Rd = q_gas / q_eps;

    if (method == Method::ReversibleAdiabat ||
        method == Method::PseudoAdiabat) {
      chi[rk] = calDlnTDlnP(*air, latent);
    } else if (method == Method::DryAdiabat) {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = calDlnTDlnP(*air, latent);
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
  EquilibrateTP(air);
  if (method != Method::ReversibleAdiabat) {
    for (int j = 0; j < NCLOUD; ++j) air->c[j] = 0;
  }
}
