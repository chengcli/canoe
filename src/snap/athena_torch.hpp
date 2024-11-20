#pragma once

// C/C++
#include <exception>

// athena
#include <athena/athena_arrays.hpp>

// torch
#include <torch/torch.h>

torch::Tensor to_torch(AthenaArray<double> const& w) {
  auto nx1 = w.GetDim1();
  auto nx2 = w.GetDim2();
  auto nx3 = w.GetDim3();
  auto nx4 = w.GetDim4();
  auto nx5 = w.GetDim5();
  auto nx6 = w.GetDim6();

  if (nx5 != 1 || nx6 != 1) {
    throw std::runtime_error("to_torch: nx5 and nx6 must be 1");
  }

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;
  int64_t str4 = nx3 * nx2 * nx1;

  return torch::from_blob(const_cast<Real*>(w.data()), {nx4, nx3, nx2, nx1},
                          {str4, str3, str2, str1}, nullptr,
                          torch::dtype(torch::kDouble));
}

AthenaArray<double> to_athena(torch::Tensor tensor) {
  auto nx1 = tensor.size(3);
  auto nx2 = tensor.size(2);
  auto nx3 = tensor.size(1);
  auto nx4 = tensor.size(0);

  int64_t str1 = 1;
  int64_t str2 = nx1;
  int64_t str3 = nx2 * nx1;
  int64_t str4 = nx3 * nx2 * nx1;

  AthenaArray<double> p(nx4, nx3, nx2, nx1);

  // create a temporary tensor holder
  torch::Tensor tmp =
      torch::from_blob(p.data(), {nx4, nx3, nx2, nx1}, {str4, str3, str2, str1},
                       nullptr, torch::dtype(torch::kDouble));
  tmp.copy_(tensor);

  return p;
}
