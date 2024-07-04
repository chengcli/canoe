#pragma once

#include <torch/torch.h>

#include <vector>

class Interpolation {
 public:
  Interpolation() {}
  virtual ~Interpolation() {}

  void ToDevice(c10::DeviceType dtype) {
    for (auto& c : cm_) c = c.to(dtype);
    for (auto& c : cp_) c = c.to(dtype);
  }

  virtual torch::Tensor left(torch::Tensor const& phi, std::string pat) = 0;
  virtual torch::Tensor right(torch::Tensor const& phi, std::string pat) = 0;

 protected:
  std::vector<torch::Tensor> cm_;
  std::vector<torch::Tensor> cp_;
};

class Center5Interp : public Interpolation {
 public:
  Center5Interp(c10::DeviceType dtype = c10::kCPU);

  torch::Tensor left(torch::Tensor const& phi, std::string pat) override;
  torch::Tensor right(torch::Tensor const& phi, std::string pat) override;
};

class Weno5Interp : public Interpolation {
 public:
  Weno5Interp(c10::DeviceType dtype = c10::kCPU);

  torch::Tensor left(torch::Tensor const& phi, std::string pat) override;
  torch::Tensor right(torch::Tensor const& phi, std::string pat) override;
};
