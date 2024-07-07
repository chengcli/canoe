#pragma once

#include <torch/torch.h>

// torch module macro is defined in
// torch/csrc/api/include/torch/nn/pimpl.h

namespace canoe {

class HydrodynamicsImpl : public torch::nn::Cloneable<HydrodynamicsImpl> {
 public:
  // Constructor to initialize the layers
  Hydro(int64_t nhydro, int64_t nghost, int64_t nx1, int64_t nx2, int64_t nx3) {
    gammad_ref_ = register_parameter("gammad_ref", torch::ones({1}) * 1.4);

    // Construct and register two Linear submodules
    fc1 = register_module("fc1", torch::nn::Linear(4, 5));
    fc2 = register_module("fc2", torch::nn::Linear(5, 3));

    auto ncells1 = nx1 > 1 ? nx1 + 2 * nghost : 1;
    auto ncells2 = nx2 > 1 ? nx2 + 2 * nghost : 1;
    auto ncells3 = nx3 > 1 ? nx3 + 2 * nghost : 1;

    register_buffer("u0", torch::empty({nhydro, nx3, nx2, nx1}));
    register_buffer("w", torch::empty({nhydro, nx3, nx2, nx1}));
    register_buffer("gammad", torch::empty({nx3, nx2, nx1}));
  }

  // Implement the one stage forward computation
  torch::Tensor forward_one_stage(torch::Tensor u, double dt) {
    auto w = eos_cons2prim_hydro_ideal(u, gammad);

    if (ncells1_ > 0.) {
      auto [wl, wr] = recon_weno5_hydro(w, IVX, 1);
      flux = rs_hydro_lmars(DIM1, wl, wr, gammad);
    }

    return add_flux_divergence(wght, flux, area, vol);

    // Use one of many tensor manipulation functions
    auto prim = torch::relu(fc1->forward(x));
    x = torch::dropout(x, /*p=*/0.5, /*train=*/is_training());
    x = fc2->forward(x);
    return torch::log_softmax(x, /*dim=*/1);
  }

 protected:
  // Use one of many "standard library" modules
  torch::nn::Linear fc1{nullptr}, fc2{nullptr};

  double gammad_ref_;
  bool is_physical_boundary_[6];
};

/// A `ModuleHolder` subclass for `HydrodynamicsImpl`.
/// See the documentation for `HydrodynamicsImpl` class to learn what methods it
/// provides, and examples of how to use `Hydrodynamics` with
/// `torch::nn::HydrodynamicsOptions`. See the documentation for `ModuleHolder`
/// to learn about PyTorch's module storage semantics.
TORCH_MODULE(Hydrodynamics);

}  // namespace canoe
