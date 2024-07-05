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

namespace eos {
namespace internal {

void _apply_conserved_limiter_inplace(torch::Tensor &cons);
void _apply_primitive_limiter_inplace(torch::Tensor &prim);

}  // namespace internal
}  // namespace eos

torch::Tensor eos_cons2prim_hydro_ideal(
    int64_t IVX, torch::Tensor &cons, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu = std::nullopt,
    std::optional<torch::Tensor> rcv = std::nullopt, int64_t ncloud = 0,
    std::optional<torch::TensorList> cos_theta = std::nullopt);

torch::Tensor eos_prim2cons_hydro_ideal(
    int64_t IVX, torch::Tensor &prim, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu = std::nullopt,
    std::optional<torch::Tensor> rcv = std::nullopt, int64_t ncloud = 0,
    std::optional<torch::TensorList> cos_theta = std::nullopt);
