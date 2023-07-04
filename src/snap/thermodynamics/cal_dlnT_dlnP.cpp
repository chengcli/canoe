// canoe
#include <variables.hpp>

// snap
#include "thermodynamics.hpp"

Real Thermodynamics::CaldlnTdlnP(Variable const& qfrac, int isat[]) const {
  // calculate gammad
  Real gammad = getGammad(qfrac);

  Real q_gas = 1.;
  for (int n = 0; n < NCLOUD; ++n) q_gas -= qfrac.c[n];

  Real f_sig = 1.;
  // vapor
  for (int n = 1; n <= NVAPOR; ++n)
    f_sig += qfrac.w[n] * (cp_ratio_mole_[n] - 1.);
  // cloud
  for (int n = 0; n < NCLOUD; ++n)
    f_sig += qfrac.c[n] * (cp_ratio_mole_[n] - 1.);
  Real cphat_ov_r = gammad / (gammad - 1.) * f_sig / q_gas;

  // vapor
  Real xd = 1.;
  for (int n = 1; n <= NVAPOR; ++n) xd -= qfrac.w[n];
  // clouds
  for (int n = 0; n < NCLOUD; ++n) xd -= qfrac.c[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    if (isat[iv] > 0) {
      Real latent = GetLatentHeat(iv, qfrac.w[IDN]);
      c1 += qfrac.w[iv] / xd * latent;
      c2 += qfrac.w[iv] / xd * latent * latent;
    }
    c3 += qfrac.w[iv] / xd;
  }

  return (1. + c1) / (cphat_ov_r + (c2 + c1 * c1) / (1. + c3));
}
