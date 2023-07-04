// canoe
#include <variables.hpp>

// snap
#include "thermodynamics.hpp"

Real Thermodynamics::TryEquilibriumTP(Variable const& qfrac, int ivapor,
                                      Real alpha, bool no_cloud) {
  int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
  int ic = NHYDRO - NVAPOR + nc - 1;

  Real xv = qfrac.w[ivapor];
  Real t = qfrac.w[IDN] / t3_[ivapor];
  auto svp_func = GetSatVaporPresFunc(ivapor, icloud);
  Real s = svp_func(qfrac) / qfrac.w[IPR];

  if (no_cloud) return xv - s;

  Real xc = qfrac.c[icloud];
  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real g = 1.;
  for (int n = 0; n < NCLOUD; ++n) g -= qfrac.c[n];
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
