// torch
#include <torch/torch.h>

namespace canoe::eos {
void _check_dim1(int64_t IVX, int64_t ncloud, const torch::Tensor &var,
                 const torch::Tensor &gammad) {
  if (IVX < 1) {
    AT_ERROR("IVX must be at least 1");
  }

  if (var.dim() != 4) {
    AT_ERROR("variable must have 4 dimensions");
  }

  if (gammad.dim() != 3) {
    AT_ERROR("gammad must have 3 dimensions");
  }

  if (ncloud >= IVX) {
    AT_ERROR("ncloud must be less than IVX");
  }
}

void _check_dim2(int64_t IVX, std::optional<torch::Tensor> var) {
  if (var.value().dim() != 1 || var.value().size(0) != IVX - 1) {
    AT_ERROR("array must have size", IVX - 1, ", and dimension 1");
  }
}

//! \todo Implement this function.
void _apply_conserved_limiter_inplace(torch::Tensor &cons) {}

//! \todo Implement this function.
void _apply_primitive_limiter_inplace(torch::Tensor &prim) {}
}  // namespace canoe::eos
