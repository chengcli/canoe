// C/C++
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// canoe
#include <configure.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/vectorize.hpp>

// thermodynamics
#include "thermodynamics_helper.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/water_vapors.hpp"

void __attribute__((weak)) update_gamma(Real *gamma, Real const q[]) {}

/*Real qhat_eps(Real const q[], Real const eps[])
{
  Real feps = 1.;
  for (int n = 1; n < NMASS; ++n)
    feps += q[n]*(eps[n] - 1.);
  return q_gas(q)/feps;
}*/

Real VaporCloudEquilibrium(Real const q[], int iv, int ic, Real t3, Real p3,
                           Real alpha, Real beta, Real delta, bool no_cloud) {
  Real xv = q[iv];
  Real t = q[IDN] / t3;
  Real s;
  if (iv == AMMONIA_VAPOR_ID)
    s = sat_vapor_p_NH3_BriggsS(q[IDN]);
  else if (iv == WATER_VAPOR_ID)
    s = sat_vapor_p_H2O_BriggsS(q[IDN]);
  else
    s = SatVaporPresIdeal(t, p3, beta, delta);
  s /= q[IPR];

  if (no_cloud) return xv - s;

  Real xc = q[ic];
  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real g = 1.;
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) g -= q[n];
  g -= xv;

  Real s1 = s / (1. - s);
  Real rate =
      (xv - s1 * g) / (1. + alpha * g * (beta / t - delta) * s1 / (1. - s));

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  // evaporate at most xc cloud
  if (rate < 0.) rate = std::min(0., std::max(rate, -xc));
  return rate;
}

Real dlnTdlnP(Real const q[], int const isat[], Real const rcp[],
              Real const beta[], Real const delta[], Real const t3[],
              Real gamma) {
  Real q_gas = 1.;
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) q_gas -= q[n];

  Real f_sig = 1.;
  // vapor
  for (int n = 1; n <= NVAPOR; ++n) f_sig += q[n] * (rcp[n] - 1.);
  // cloud
  for (int n = 0; n < 2 * NVAPOR; ++n)
    f_sig += q[NHYDRO + n] * (rcp[1 + NVAPOR + n] - 1.);
  Real cphat_ov_r = gamma / (gamma - 1.) * f_sig / q_gas;

  // vapor
  Real xd = 1.;
  for (int n = 1; n <= NVAPOR; ++n) xd -= q[n];
  // clouds
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) xd -= q[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    if (isat[iv] > 0) {
      int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
      // int nc = iv + NVAPOR;
      Real latent = beta[nc] * t3[iv] / q[IDN] - delta[nc];
      c1 += q[iv] / xd * latent;
      c2 += q[iv] / xd * latent * latent;
    }
    c3 += q[iv] / xd;
  }

  return (1. + c1) / (cphat_ov_r + (c2 + c1 * c1) / (1. + c3));
}
