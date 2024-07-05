#pragma once

namespace c10 {
template <typename T>
class ArrayRef;
}

namespace at {
class Tensor;
}

namespace torch {
using Tensor = at::Tensor;
using TensorList = c10::ArrayRef<Tensor>;
}  // namespace torch

void add_flux_divergence_inplace(double wght, torch::TensorList flux,
                                 torch::TensorList area,
                                 torch::Tensor const& vol, torch::Tensor& out);
