#pragma once

namespace at {
class Tensor;
}  // namespace at

namespace torch {
using Tensor = at::Tensor;
}  // namespace torch

namespace canoe {
std::pair<torch::Tensor, torch::Tensor, torch::Tensor> decom_obtain_anomaly(
    int64_t is, int64_t ie, torch::Tensor const& w, torch::Tensor const& dx1f,
    torch::Tensor const& grav);

void decom_apply_anomaly_inplace(torch::Tensor const& psf,
                                 torch::Tensor const& tsf, torch::Tensor& wl,
                                 torch::Tensor& wr)
}  // namespace canoe
