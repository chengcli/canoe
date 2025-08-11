#pragma once

// athena
#include <athena/athena.hpp>

// snap
#include <snap/snap.h>

#include <snap/eos/ideal_moist.hpp>

// canoe
#include "hydro.hpp"

inline static Real get_rd() {
  return kintera::constants::Rgas / kintera::species_weights[0];
}

inline static Real get_gammad() {
  auto cvd = kintera::species_cref_R[0];
  return (cvd + 1.) / cvd;
}

inline torch::Tensor get_temp(snap::IdealMoist peos,
                              AthenaArray<Real> const& hydro_w) {
  auto w = get_all(hydro_w);
  return peos->compute("W->T", {w});
}

inline torch::Tensor get_cv(snap::IdealMoist peos,
                            AthenaArray<Real> const& hydro_w) {
  auto pthermo = peos->pthermo;
  auto mud = kintera::species_weights[0];
  auto Rd = kintera::constants::Rgas / mud;
  auto cvd = kintera::species_cref_R[0] * Rd;
  int ny = pthermo->options.vapor_ids().size() +
           pthermo->options.cloud_ids().size() - 1;
  auto w = get_all(hydro_w);
  auto fsig = peos->f_sig(w.narrow(0, snap::ICY, ny));
  return cvd * fsig;
}

//! \brief Effective polytropic index
//!
//! Eq.63 in Li2019
//! $\gamma = \frac{c_p}{c_v}$
//! \return $\gamma$
template <typename T>
Real get_gamma(snap::IdealMoist const& peos, T u, Real rho = 1.) {
  Real fsig = 1., feps = 1.;
  auto cv_ratio_m1 = peos->cv_ratio_m1.accessor<Real, 1>();
  auto inv_mu_ratio_m1 = peos->inv_mu_ratio_m1.accessor<Real, 1>();

  int nvapor = peos->pthermo->options.vapor_ids().size() - 1;
  int ncloud = peos->pthermo->options.cloud_ids().size();

  for (int n = 1; n <= nvapor; ++n) {
    fsig += u[n] / rho * cv_ratio_m1[n - 1];
    feps += u[n] / rho * inv_mu_ratio_m1[n - 1];
  }

  for (int n = 1 + nvapor; n < 1 + nvapor + ncloud; ++n) {
    fsig += u[n] / rho * cv_ratio_m1[n - 1];
    feps -= u[n] / rho;
  }

  return 1. + (get_gammad() - 1.) * feps / fsig;
}
