// C/C++
#include <algorithm>
#include <iomanip>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

// TODO(cli) : check
void Thermodynamics::rk4IntegrateZ(Variable *qfrac, Real dz, Method method,
                                   Real grav, Real adTdz) const {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = qfrac->w[IDN];
  Real pres = qfrac->w[IPR];
  Real dTdz[4], chi[4];
  Real latent[1 + NVAPOR];

  for (int rk = 0; rk < 4; ++rk) {
    EquilibrateTP(qfrac);
    if (method != Method::ReversibleAdiabat)
      for (int j = 0; j < NCLOUD; ++j) qfrac->c[j] = 0;

    for (int i = 1; i <= NVAPOR; ++i) {
      // make a trial run to get latent heat
      qfrac->w[i] += 1.E-6;
      auto rates = TryEquilibriumTP_VaporCloud(*qfrac, i);
      latent[i] = GetLatentHeatMole(i, rates, qfrac->w[IDN]) /
                  (Constants::Rgas * qfrac->w[IDN]);
      qfrac->w[i] -= 1.E-6;
    }

    Real q_gas = 1., q_eps = 1.;
    for (int i = 1; i <= NVAPOR; ++i) {
      q_eps += qfrac->w[i] * (mu_ratio_[i] - 1.);
    }

    for (int j = 0; j < NCLOUD; ++j) {
      q_eps += qfrac->c[j] * (mu_ratio_[1 + NVAPOR + j] - 1.);
      q_gas += -qfrac->c[j];
    }

    Real g_ov_Rd = grav / Rd_;
    Real R_ov_Rd = q_gas / q_eps;

    if (method == Method::ReversibleAdiabat ||
        method == Method::PseudoAdiabat) {
      chi[rk] = calDlnTDlnP(*qfrac, latent);
    } else if (method == Method::DryAdiabat) {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = calDlnTDlnP(*qfrac, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }

    dTdz[rk] = -chi[rk] * g_ov_Rd / R_ov_Rd + adTdz;
    chi[rk] = -R_ov_Rd / g_ov_Rd * dTdz[rk];

    // integrate over dz
    Real chi_avg;
    if (rk < 3) {
      qfrac->w[IDN] = temp + dTdz[rk] * dz * step[rk];
      chi_avg = chi[rk];
    } else {
      qfrac->w[IDN] =
          temp +
          1. / 6. * (dTdz[0] + 2. * dTdz[1] + 2. * dTdz[2] + dTdz[3]) * dz;
      chi_avg = 1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]);
    }

    if (!(qfrac->w[IDN] > 0.)) qfrac->w[IDN] = temp;

    if (fabs(qfrac->w[IDN] - temp) > 0.01) {
      qfrac->w[IPR] = pres * pow(qfrac->w[IDN] / temp, 1. / chi_avg);
    } else {  // isothermal limit
      qfrac->w[IPR] =
          pres * exp(-2. * g_ov_Rd * dz / (R_ov_Rd * (qfrac->w[IDN] + temp)));
    }
  }

  // recondensation
  EquilibrateTP(qfrac);
  if (method != Method::ReversibleAdiabat)
    for (int j = 0; j < NCLOUD; ++j) qfrac->c[j] = 0;
}
