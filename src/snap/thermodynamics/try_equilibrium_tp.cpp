// C/C++
#include <algorithm>

// canoe
#include <variable.hpp>

// snap
#include "thermodynamics.hpp"

std::vector<Real> Thermodynamics::TryEquilibriumTP(Variable const& qfrac,
                                                   int ivapor, Real l_over_cv,
                                                   bool misty) const {
  Real xv = qfrac.w[ivapor];
  Real t = qfrac.w[IDN] / t3_[ivapor];
  std::vector<Real> rates(1 + cloud_index_set_[ivapor].size(), 0.);

  for (int n = 0; n < cloud_index_set_[ivapor].size(); ++n) {
    int jcloud = cloud_index_set_[ivapor][n];
    Real xs = svp_func_[ivapor][n](qfrac, ivapor, jcloud) / qfrac.w[IPR];

    Real xc = qfrac.c[jcloud];

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
    Real rate = (s1 * g - xv) /
                (1. + l_over_cv * g * (beta_[jcloud] / t - delta_[jcloud]) *
                          s1 / (1. - xs));

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
