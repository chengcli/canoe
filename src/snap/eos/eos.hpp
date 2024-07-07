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

namespace canoe {
namespace eos {
void _check_dim1(int64_t IVX, int64_t ncloud, const torch::Tensor &var,
                 const torch::Tensor &gammad);
void _check_dim2(int64_t IVX, std::optional<torch::Tensor> var);

void _apply_conserved_limiter_inplace(torch::Tensor &cons);
void _apply_primitive_limiter_inplace(torch::Tensor &prim);
}  // namespace eos

torch::Tensor eos_cons2prim_hydro_ideal(
    torch::Tensor &cons, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu = std::nullopt,
    std::optional<torch::Tensor> rcv = std::nullopt, int64_t ncloud = 0,
    std::optional<torch::TensorList> cos_theta = std::nullopt);

torch::Tensor eos_prim2cons_hydro_ideal(
    torch::Tensor &prim, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu = std::nullopt,
    std::optional<torch::Tensor> rcv = std::nullopt, int64_t ncloud = 0,
    std::optional<torch::TensorList> cos_theta = std::nullopt);
}  // namespace canoe
