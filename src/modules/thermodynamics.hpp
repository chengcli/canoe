#pragma once

#include <torch/arg.h>
#include <torch/nn.h>

namespace canoe {

class ThermodynamicsOptions {
 public:
  ThermodynamicsOptions() {}

  TORCH_ARG(double, gammad_ref) = 1.4;
  TORCH_ARG(int64_t, nvapor) = 0;
  TORCH_ARG(int64_t, ncloud) = 0;
};

class ThermodynamicsImpl : public torch::nn::Cloneable<ThermodynamicsImpl> {
 public:
  //! options with which this `Thermodynamics` was constructed
  ThermodynamicsOptions options;

  // Constructor to initialize the layers
  explicit ThermodynamicsImpl(const ThermodynamicsOptions& options_);

  void reset() override;

  int64_t nvapor() const;
  int64_t ncloud() const;

  torch::Tensor mu_ratio() const;
  torch::Tensor cv_ratio_mass() const;
  torch::Tensor gammad(torch::Tensor const& hydro_u) const;

  // Implement the one stage forward computation
  std::pair<torch::Tensor, torch::Tensor> forward(torch::Tensor hydro_u,
                                                  torch::Tensor scalar_u);
};

TORCH_MODULE(Thermodynamics);

}  // namespace canoe
