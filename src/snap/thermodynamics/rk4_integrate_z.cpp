// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

// TODO(cli) : check
void Thermodynamics::rk4IntegrateZ(Variable *qfrac, Real dz, Method method,
                                   Real g_ov_Rd, Real adTdz) const {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = qfrac->w[IDN];
  Real pres = qfrac->w[IPR];
  Real dTdz[4], chi[4];
  Real latent[1 + NVAPOR];

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor and cloud
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      setTotalEquivalentVapor(qfrac, iv);
      latent[iv] = 0;
    }

    Real q_gas = 1., q_eps = 1.;
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      auto rates = TryEquilibriumTP(*qfrac, iv);

      // saturation indicator
      latent[iv] = GetLatentHeatMole(iv, rates, qfrac->w[IDN]);

      // vapor condensation rate
      qfrac->w[iv] += rates[0];
      q_eps += rates[0] * (mu_ratio_[iv] - 1.);

      // cloud concentration rates
      if (method == Method::ReversibleAdiabat) {
        for (int n = 1; n < rates.size(); ++n) {
          qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
          q_eps += rates[n] * (mu_ratio_[iv] - 1.);
          q_gas -= rates[n];
        }
      }
    }
    Real R_ov_Rd = q_gas / q_eps;

    if (method == Method::ReversibleAdiabat ||
        method == Method::PseudoAdiabat) {
      chi[rk] = calDlnTDlnP(*qfrac, latent);
    } else if (method == Method::DryAdiabat) {
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
    if (fabs(qfrac->w[IDN] - temp) > 0.1)  // isothermal limit
      qfrac->w[IPR] = pres * pow(qfrac->w[IDN] / temp, 1. / chi_avg);
    else
      qfrac->w[IPR] =
          pres * exp(-2. * g_ov_Rd * dz / (R_ov_Rd * (qfrac->w[IDN] + temp)));
  }

  // recondensation
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    setTotalEquivalentVapor(qfrac, iv);
    auto rates = TryEquilibriumTP(*qfrac, iv);

    // vapor condensation rate
    qfrac->w[iv] += rates[0];

    // cloud concentration rates
    if (method == Method::ReversibleAdiabat) {
      for (int n = 1; n < rates.size(); ++n)
        qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
    }
  }
}
