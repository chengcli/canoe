#pragma once

// athena
#include <athena/athena.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// canoe
#include "hydro.hpp"

inline AthenaArray<Real> get_temp(snap::IdealMoist peos,
                                  AthenaArray<Real> const& hydro_w) {
  auto w = get_all(hydro_w);
  auto T = peos->compute("W->T", {w});

  auto nx1 = T.size(2);
  auto nx2 = T.size(1);
  auto nx3 = T.size(0);

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;

  AthenaArray<double> p(nx3, nx2, nx1);

  // create a temporary tensor holder
  torch::Tensor tmp =
      torch::from_blob(p.data(), {nx3, nx2, nx1}, {str3, str2, str1}, nullptr,
                       torch::dtype(torch::kDouble));
  tmp.copy_(T);

  return p;
}
