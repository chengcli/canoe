#include "thermodynamics.hpp"

//! Eq.1 in Li2018
Real cal_dlnT_dlnP(Real const* xfrac, Real gammad, Real const* cp_ratio_mole,
                   Real const* latent) {
  Real q_gas = 1.;
  for (int n = 1 + NVAPOR; n <= NVAPOR + NCLOUD; ++n) q_gas -= xfrac[n];

  Real f_sig = 1.;
  // vapor
  for (int n = 1; n <= NVAPOR + NCLOUD; ++n) {
    f_sig += xfrac[n] * (cp_ratio_mole[n] - 1.);
  }
  Real cphat_ov_r = gammad / (gammad - 1.) * f_sig / q_gas;

  // vapor
  Real xd = q_gas;
  for (int n = 1; n <= NVAPOR; ++n) xd -= xfrac[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    c1 += xfrac[iv] / xd * latent[iv];
    c2 += xfrac[iv] / xd * latent[iv] * latent[iv];
    c3 += xfrac[iv] / xd;
  }

  return (1. + c1) / (cphat_ov_r + (c2 + c1 * c1) / (1. + c3));
}
