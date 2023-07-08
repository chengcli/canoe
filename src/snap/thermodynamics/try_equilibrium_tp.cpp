// C/C++
#include <algorithm>

// canoe
#include <variable.hpp>

// snap
#include "thermodynamics.hpp"

std::vector<Real> Thermodynamics::TryEquilibriumTP(Variable const& qfrac, int i,
                                                   Real l_over_cv,
                                                   bool misty) const {
  Real xv = qfrac.w[i];
  Real t = qfrac.w[IDN] / t3_[i];
  std::vector<Real> rates(1 + cloud_index_set_[i].size(), 0.);

  for (int n = 0; n < cloud_index_set_[i].size(); ++n) {
    int j = cloud_index_set_[i][n];
    Real xs = svp_func_[i][n](qfrac, i, j) / qfrac.w[IPR];
    Real xc = qfrac.c[j];

    if (misty) {  // no cloud in variable
      rates[0] += xs - xv;
      continue;
    }

    // if saturation vapor pressure is larger than the total pressure
    // evaporate all condensates
    if (xs > 1.) {
      rates[0] += xc;
      rates[1 + n] = -xc;
      continue;
    }

    Real g = 1.;
    for (int n = 0; n < NCLOUD; ++n) g -= qfrac.c[n];
    g -= xv;

    Real s1 = xs / (1. - xs);
    Real rate =
        (s1 * g - xv) /
        (1. + l_over_cv * g * (beta_[j] / t - delta_[j]) * s1 / (1. - xs));

    // condensate at most xv vapor
    if (rate < 0.) {
      rates[0] += -std::min(std::abs(rate), xv);
      rates[1 + n] = std::min(std::abs(rate), xv);
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
