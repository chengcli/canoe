// C/C++
#include <map>
#include <string>

// torch
#include <ATen/Tensor.h>

// harp
#include <harp/radiation/radiation.hpp>

struct RadiationData {
  double counter = 0.0;
  double cooldown = 0.0;
  at::Tensor net_flux;
  std::map<std::string, at::Tensor> rad_bc;
};
