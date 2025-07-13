#pragma once

// torch
#include <torch/torch.h>

// athena
#include <athena/athena.hpp>

torch::Tensor get_dens(AthenaArray<Real> const& w) {
  int64_t str1 = 1;
  int64_t str2 = w.GetDim1();
  int64_t str3 = w.GetDim2() * w.GetDim1();
  int64_t str4 = w.GetDim3() * w.GetDim2() * w.GetDim3();

  int64_t n1 = w.GetDim1();
  int64_t n2 = w.GetDim2();
  int64_t n3 = w.GetDim3();

  return torch::from_blob(w.data(), {n3, n2, n1}, {str3, str2, str1}, nullptr,
                          torch::dtype(torch::kFloat64));
}

torch::Tensor get_pres(AthenaArray<Real> const& w) {
  int64_t str1 = 1;
  int64_t str2 = w.GetDim1();
  int64_t str3 = w.GetDim2() * w.GetDim1();
  int64_t str4 = w.GetDim3() * w.GetDim2() * w.GetDim3();

  int64_t n1 = w.GetDim1();
  int64_t n2 = w.GetDim2();
  int64_t n3 = w.GetDim3();

  return torch::from_blob(w.data() + IPR * str4, {n3, n2, n1},
                          {str3, str2, str1}, nullptr,
                          torch::dtype(torch::kFloat64));
}

torch::Tensor get_yfrac(AthenaArray<Real> const& w) {
  int64_t str1 = 1;
  int64_t str2 = w.GetDim1();
  int64_t str3 = w.GetDim2() * w.GetDim1();
  int64_t str4 = w.GetDim3() * w.GetDim2() * w.GetDim3();

  int64_t n1 = w.GetDim1();
  int64_t n2 = w.GetDim2();
  int64_t n3 = w.GetDim3();
  int64_t n4 = IVX - 1;

  return torch::from_blob(w.data() + str4, {n4, n3, n2, n1},
                          {str4, str3, str2, str1}, nullptr,
                          torch::dtype(torch::kFloat64));
}

Real& get_val(torch::Tensor w, int k, int j, int i) {
  auto a = w.accessor<Real, 3>();
  return a[k][j][i];
}

Real& get_val(torch::Tensor w, int n, int k, int j, int i) {
  auto a = w.accessor<Real, 4>();
  return a[n][k][j][i];
}
