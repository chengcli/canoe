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

void vec_lower_inplace(torch::Tensor &prim, torch::TensorList cos_theta);
void vec_raise_inplace(torch::Tensor &prim, torch::TensorList cos_theta);
