// torch
#include "hydrodynamics.hpp"

#include <torch/torch.h>

#include <torch/eos/eos.hpp>
#include <torch/hydro/hydro.hpp>
#include <torch/recon/recon.hpp>
#include <torch/riemann/riemann.hpp>

enum { DIM1 = 3, DIM2 = 2, DIM3 = 1 };

namespace canoe {

HydrodynamicsImpl::HydrodynamicsImpl(const HydrodynamicsOptions& options_,
                                     Coordinates coord,
                                     std::optional<EquationOfState> eos) {
  coord_ = register_module("coord", coord);

  if (eos.has_value()) {
    eos_ = register_module("eos", eos.value());
  } else {
    eos_ = register_module("eos", EquationOfState(EquationOfStateOptions()));
  }

  reset();
}

void HydrodynamicsImpl::reset() {
  register_buffer(
      "flux1", torch::zeros({eos_->NHYDRO, coord_->ncells3(), coord_->ncells2(),
                             coord_->ncells1() + 1}));
  register_buffer("flux2",
                  torch::zeros({eos_->NHYDRO, coord_->ncells3(),
                                coord_->ncells2() + 1, coord_->ncells1()}));
  register_buffer("flux3",
                  torch::zeros({eos_->NHYDRO, coord_->ncells3() + 1,
                                coord_->ncells2(), coord_->ncells1()}));
  register_buffer("w", torch::zeros({eos_->NHYDRO, coord_->ncells3(),
                                     coord_->ncells2(), coord_->ncells1()}));
}

float HydrodynamicsImpl::cfl() const { return options.cfl(); }

torch::Tensor HydrodynamicsImpl::forward(torch::Tensor u, double dt) {
  auto buf = named_buffers();
  auto gammad = eos_->gammad(u);

  buf["w"] = eos_cons2prim_hydro_ideal(u, gammad, eos_->rmu(), eos_->rcv(),
                                       eos_->ncloud(), coord_->cos_theta());

  // bc_apply_boundary_conditions_inplace(buf["w"], coord_.options());

  if (u.size(DIM1) > 0.) {
    auto result = recon_weno5_hydro(buf["w"], eos_->IVX, 1);
    buf["flux1"] =
        rs_hydro_lmars(DIM1, result.first, result.first, gammad, eos_->rmu(),
                       eos_->rcv(), eos_->ncloud(), coord_->cos_theta());
  }

  // auto dt = max_timestep(w);

  return u - cfl() * dt *
                 flux_divergence({buf["flux1"], buf["flux2"], buf["flux3"]},
                                 coord_->areas(), coord_->vol());

  // Use one of many tensor manipulation functions
  // auto prim = torch::relu(fc1->forward(x));
  // x = torch::dropout(x, /*p=*/0.5, /*train=*/is_training());
  // x = fc2->forward(x);
  // return torch::log_softmax(x, /*dim=*/1);
}

}  // namespace canoe
