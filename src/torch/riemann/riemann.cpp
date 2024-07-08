// torch
#include <torch/torch.h>

// snap
#include "riemann.hpp"

namespace canoe::riemann {

//! \todo Implement this function.
void _prim2local_inplace(int64_t ivx, torch::Tensor const& w,
                         torch::TensorList cos_theta) {}

//! \todo Implement this function.
void _flux2global_inplace(int64_t ivx, torch::Tensor const& w,
                          torch::TensorList cos_theta) {}

}  // namespace canoe::riemann
