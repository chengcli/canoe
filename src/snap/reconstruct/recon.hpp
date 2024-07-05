// torch
#include <torch/torch.h>

std::pair<torch::Tensor, torch::Tensor> recon_weno5_hydro(
    const torch::Tensor &w, int64_t IVX, int64_t dim, bool mixed = true,
    bool is_boundary_lower = false, bool is_boundary_upper = false);
