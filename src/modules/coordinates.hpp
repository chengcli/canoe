#pragma once

#include <torch/arg.h>
#include <torch/nn.h>

namespace canoe {

struct CoordinatesOptions {
  CoordinatesOptions() {}

  TORCH_ARG(int64_t, nx1) = 1;
  TORCH_ARG(int64_t, nx2) = 1;
  TORCH_ARG(int64_t, nx3) = 1;
  TORCH_ARG(int64_t, nghost) = 0;
};

class CoordinatesImpl : public torch::nn::Cloneable<CoordinatesImpl> {
 public:
  //! options with which this `Coordinates` was constructed
  CoordinatesOptions options;

  CoordinatesImpl();

  // Constructor to initialize the layers
  explicit CoordinatesImpl(const CoordinatesOptions& options_);

  void reset() override;

  int64_t ncells3() const;
  int64_t ncells2() const;
  int64_t ncells1() const;

  torch::Tensor vol() const;
  torch::TensorList areas() const;
  torch::Tensor areas(int64_t dim) const;

  torch::TensorList cos_theta() const;
  torch::Tensor cos_theta(int64_t dim) const;

  bool is_physical_boundary_low(int64_t dim) const;
  bool is_physical_boundary_high(int64_t dim) const;

 protected:
  bool is_physical_boundary_[6];
};

TORCH_MODULE(Coordinates);

}  // namespace canoe
