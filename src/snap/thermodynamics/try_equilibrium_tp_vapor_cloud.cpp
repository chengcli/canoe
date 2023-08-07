// C/C++
#include <algorithm>

// canoe
#include <variable.hpp>

// snap
#include "thermodynamics.hpp"

// Calculates phase equilibrium of
// Vapor <=> Cloud
//
// Example phase equilibrium:
// H2O -> H2O(l)
//
RealArrayX Thermodynamics::TryEquilibriumTP_VaporCloud(Variable const& qfrac,
                                                       int i, Real cv_hat,
                                                       bool misty) const {
  Real xv = qfrac.w[i], xg = 1.;
  Real t = qfrac.w[IDN] / t3_[i];
  std::vector<Real> rates(1 + cloud_index_set_[i].size(), 0.);

#pragma omp simd reduction(+ : xg)
  for (int n = 0; n < NCLOUD; ++n) xg += -qfrac.c[n];

  for (int n = 0; n < cloud_index_set_[i].size(); ++n) {
    int j = cloud_index_set_[i][n];
    Real xs = svp_func1_[i][n](qfrac, i, j) / qfrac.w[IPR];
    Real xc = qfrac.c[j];

    if (misty) {  // in a cloudy ambient environment
      rates[0] += xs - xv / xg;
      continue;
    }

    // if saturation vapor pressure is larger than the total pressure
    // evaporate all condensates
    if (xs > 1.) {
      rates[0] += xc;
      rates[1 + n] = -xc;
      continue;
    }

    Real alpha = 0.;
    xg -= xv;

    Real lv = beta_[1 + NVAPOR + j] / t - delta_[1 + NVAPOR + j];
    if (cv_hat > 0.) alpha = (lv - 1.) / cv_hat;

    Real s1 = xs / (1. - xs);
    Real rate = (s1 * xg - xv) / (1. + alpha * xg * lv * s1 / (1. - xs));

    // condensate at most xv vapor
    if (rate < 0.) {
      rates[0] += -std::min(-rate, xv);
      rates[1 + n] = std::min(-rate, xv);
    }

    // evaporate at most xc cloud
    if (rate > 0.) {
      rates[0] += std::min(rate, xc);
      rates[1 + n] = -std::min(rate, xc);
    }
  }

  // scale total rate
  if (rates[0] < 0. && std::abs(rates[0]) > xv) {
    Real r = xv / std::abs(rates[0]);
    for (auto& rate : rates) rate *= r;
  }

  return rates;
}
