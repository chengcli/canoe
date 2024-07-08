#include "equation_of_state.hpp"

#include <torch/torch.h>

#include "thermodynamics.hpp"

namespace canoe {

EquationOfStateImpl::EquationOfStateImpl() : options(EquationOfStateOptions()) {
  reset();
}

EquationOfStateImpl::EquationOfStateImpl(const EquationOfStateOptions& options_,
                                         std::optional<Thermodynamics> thermo)
    : options(options_) {
  if (thermo.has_value()) {
    thermo_ = register_module("thermo", thermo.value());
  } else {
    thermo_ = register_module("thermo", Thermodynamics());
  }

  reset();
}

void EquationOfStateImpl::reset() {
  IVX = 1 + thermo_->nvapor() + thermo_->ncloud();
  IPR = 3 + IVX;
  NHYDRO = 1 + IPR;

  register_buffer("rmu", torch::zeros({IVX - 1}));
  register_buffer("rcv", torch::zeros({IVX - 1}));
  register_buffer("gammad", torch::zeros({IVX - 1}));
}

void EquationOfStateImpl::update_thermo(torch::Tensor const& hydro_u) const {
  auto buf = named_buffers(/*recurse=*/false);

  buf["rcv"] = (thermo_->cv_ratio_mass() - 1.).slice(0, 1, IVX);
  buf["mu_ratio"] = (1. / thermo_->mu_ratio()).slice(0, 1, IVX);
  buf["gammad"] = thermo_->gammad(hydro_u);
}

torch::Tensor EquationOfStateImpl::gammad(torch::Tensor const& hydro_u) const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["gammad"];
}

torch::Tensor EquationOfStateImpl::rmu() const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["rmu"];
}

torch::Tensor EquationOfStateImpl::rcv() const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["rcv"];
}

int64_t EquationOfStateImpl::ncloud() const { return thermo_->ncloud(); }

}  // namespace canoe
