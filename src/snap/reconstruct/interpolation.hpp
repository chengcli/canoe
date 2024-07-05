#pragma once

// C/C++
#include <vector>

// torch
#include <c10/core/DeviceType.h>

namespace at {
class Tensor;
}

namespace torch {
using Tensor = at::Tensor;
}

class Interpolation {
 public:
  Interpolation() {}
  virtual ~Interpolation() {}

  void ToDevice(c10::DeviceType dtype);
  virtual torch::Tensor left(torch::Tensor const& phi) const = 0;
  virtual torch::Tensor right(torch::Tensor const& phi) const = 0;

 protected:
  std::vector<torch::Tensor> cm_;
  std::vector<torch::Tensor> cp_;
};

class Center5Interp : public Interpolation {
 public:
  explicit Center5Interp(c10::DeviceType dtype = c10::kCPU);

  torch::Tensor left(torch::Tensor const& phi) const override;
  torch::Tensor right(torch::Tensor const& phi) const override;
};

class Weno5Interp : public Interpolation {
 public:
  explicit Weno5Interp(c10::DeviceType dtype = c10::kCPU);

  torch::Tensor left(torch::Tensor const& phi) const override;
  torch::Tensor right(torch::Tensor const& phi) const override;
};
