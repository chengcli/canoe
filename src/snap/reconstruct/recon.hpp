#pragma once

namespace at {
class Tensor;
}  // namespace at

namespace torch {
using Tensor = at::Tensor;
}  // namespace torch

namespace canoe {
std::pair<torch::Tensor, torch::Tensor> recon_weno5_hydro(
    const torch::Tensor &w, int64_t IVX, int64_t dim, bool mixed = true,
    bool is_boundary_lower = false, bool is_boundary_upper = false);

std::pair<torch::Tensor, torch::Tensor> recon_weno5_scalar(
    const torch::Tensor &w, int64_t dim);
}  // namespace canoe
