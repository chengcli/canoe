#pragma once

#include <torch/arg.h>
#include <torch/nn.h>

namespace canoe {

class Hydrodynamics;
class Transport;

class SSPRKOptions {
 public:
  SSPRKOptions();

  TORCH_ARG(int64_t, nstage) = 1;
};

class SSPRKImpl : public torch::nn::Cloneable<SSPRKImpl> {
 public:
  //! options with which this `SSPRK` was constructed
  SSPRKOptions options;

  // Constructor to initialize the layers
  explicit SSPRKImpl(const SSPRKOptions& options_);

  void reset() override;

  // Implement the one stage forward computation
  std::pair<torch::Tensor, torch::Tensor> forward(torch::Tensor hydro_u,
                                                  torch::Tensor scalar_u);

 protected:
  Hydrodynamics hydro_;
  Transport transport_;
};

}  // namespace canoe
