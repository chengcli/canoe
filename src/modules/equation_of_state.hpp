#pragma once

#include <torch/arg.h>
#include <torch/nn.h>

#include "thermodynamics.hpp"

namespace Cantera {
class Kinetics;
}  // namespace Cantera

namespace canoe {

class Thermodynamics;

struct EquationOfStateOptions {
  EquationOfStateOptions() {}
};

class EquationOfStateImpl : public torch::nn::Cloneable<EquationOfStateImpl> {
 public:
  //! options with which this `EquationOfState` was constructed
  EquationOfStateOptions options;

  int64_t IDN = 0;
  int64_t IVX = 1;
  int64_t IPR = 4;
  int64_t NHYDRO = 5;

  // Constructor to initialize the layers
  explicit EquationOfStateImpl(
      const EquationOfStateOptions& options_,
      std::optional<Thermodynamics> thermo = std::nullopt);

  void reset() override;

  void update_thermo(torch::Tensor const& hydro_u) const;

  torch::Tensor gammad(torch::Tensor const& hydro_u) const;
  torch::Tensor rmu() const;
  torch::Tensor rcv() const;
  int64_t ncloud() const;

 protected:
  Thermodynamics thermo_ = nullptr;
  std::shared_ptr<Cantera::Kinetics> kinetics_;
};

TORCH_MODULE(EquationOfState);

}  // namespace canoe
