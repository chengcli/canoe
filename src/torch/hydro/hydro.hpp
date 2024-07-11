#pragma once

namespace c10 {
template <typename T>
class ArrayRef;
}  // namespace c10

namespace at {
class Tensor;
}  // namespace at

namespace torch {
using Tensor = at::Tensor;
using TensorList = c10::ArrayRef<Tensor>;
}  // namespace torch

namespace canoe {
torch::Tensor flux_divergence(torch::TensorList flux, torch::TensorList area,
                              torch::Tensor const& vol);
}  // namespace canoe
