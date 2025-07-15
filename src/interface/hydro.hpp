#pragma once

// athena
#include <athena/athena.hpp>

// torch
#include <torch/torch.h>

inline torch::Tensor get_dens(AthenaArray<Real> const& w) {
  auto nx1 = w.GetDim1();
  auto nx2 = w.GetDim2();
  auto nx3 = w.GetDim3();

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;

  return torch::from_blob(const_cast<Real*>(w.data()), {nx3, nx2, nx1},
                          {str3, str2, str1}, nullptr,
                          torch::dtype(torch::kDouble));
}

inline torch::Tensor get_pres(AthenaArray<Real> const& w) {
  auto nx1 = w.GetDim1();
  auto nx2 = w.GetDim2();
  auto nx3 = w.GetDim3();

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;
  int64_t str4 = nx3 * nx2 * nx1;

  return torch::from_blob(const_cast<Real*>(w.data() + IPR * str4),
                          {nx3, nx2, nx1}, {str3, str2, str1}, nullptr,
                          torch::dtype(torch::kDouble));
}

inline torch::Tensor get_yfrac(AthenaArray<Real> const& w) {
  auto nx1 = w.GetDim1();
  auto nx2 = w.GetDim2();
  auto nx3 = w.GetDim3();
  auto nx4 = w.GetDim4();

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;
  int64_t str4 = nx3 * nx2 * nx1;

  return torch::from_blob(const_cast<Real*>(w.data() + str4),
                          {IVX - 1, nx3, nx2, nx1}, {str4, str3, str2, str1},
                          nullptr, torch::dtype(torch::kDouble));
}

inline torch::Tensor get_all(AthenaArray<Real> const& w) {
  auto nx1 = w.GetDim1();
  auto nx2 = w.GetDim2();
  auto nx3 = w.GetDim3();
  auto nx4 = w.GetDim4();

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;
  int64_t str4 = nx3 * nx2 * nx1;

  return torch::from_blob(const_cast<Real*>(w.data()), {nx4, nx3, nx2, nx1},
                          {str4, str3, str2, str1}, nullptr,
                          torch::dtype(torch::kDouble));
}

inline Real& get_val(torch::Tensor w, int k, int j, int i) {
  auto a = w.accessor<Real, 3>();
  return a[k][j][i];
}

inline Real& get_val(torch::Tensor w, int n, int k, int j, int i) {
  auto a = w.accessor<Real, 4>();
  return a[n][k][j][i];
}
