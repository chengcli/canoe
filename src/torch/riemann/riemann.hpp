#pragma once

// C/C++
#include <optional>

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
namespace riemann {
void _prim2local_inplace(int64_t ivx, torch::Tensor const &w,
                         torch::TensorList cos_theta);
void _flux2global_inplace(int64_t ivx, torch::Tensor const &w,
                          torch::TensorList cos_theta);
}  // namespace riemann

torch::Tensor rs_hydro_lmars(
    int64_t dim, const torch::Tensor &wl, const torch::Tensor &wr,
    const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu = std::nullopt,
    std::optional<torch::Tensor> rcv = std::nullopt, int64_t ncloud = 0,
    std::optional<torch::TensorList> cos_theta = std::nullopt);
}  // namespace canoe
