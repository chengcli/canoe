#pragma once

#include <torch/arg.h>
#include <torch/nn.h>

#include "coordinates.hpp"
#include "equation_of_state.hpp"

// torch module macro is defined in
// torch/csrc/api/include/torch/nn/pimpl.h

namespace canoe {

class EquationOfState;
class Coordinates;

struct HydrodynamicsOptions {
  HydrodynamicsOptions() {}

  TORCH_ARG(float, cfl) = 0.9;
  TORCH_ARG(float, grav1) = 0.;
};

class HydrodynamicsImpl : public torch::nn::Cloneable<HydrodynamicsImpl> {
 public:
  //! options with which this `Hydrodynamics` was constructed
  HydrodynamicsOptions options;

  // Constructor to initialize the layers
  explicit HydrodynamicsImpl(const HydrodynamicsOptions& options_,
                             Coordinates coord,
                             std::optional<EquationOfState> eos = std::nullopt);

  void reset() override;

  torch::Tensor initialize() const;

  float cfl() const;

  float max_timestep(torch::Tensor w) const;

  // Implement the one stage forward computation
  torch::Tensor forward(torch::Tensor u, double dt);

 protected:
  EquationOfState eos_ = nullptr;
  Coordinates coord_ = nullptr;

  // Use one of many "standard library" modules
  // torch::nn::Linear fc1{nullptr}, fc2{nullptr};
};

/// A `ModuleHolder` subclass for `HydrodynamicsImpl`.
/// See the documentation for `HydrodynamicsImpl` class to learn what methods it
/// provides, and examples of how to use `Hydrodynamics` with
/// `torch::nn::HydrodynamicsOptions`. See the documentation for `ModuleHolder`
/// to learn about PyTorch's module storage semantics.
TORCH_MODULE(Hydrodynamics);

}  // namespace canoe
