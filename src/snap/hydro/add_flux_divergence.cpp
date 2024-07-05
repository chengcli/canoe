// torch
#include <torch/torch.h>

// snap
#include "hydro.hpp"

enum {
  DIM1 = 3,
  DIM2 = 2,
  DIM3 = 1,
};

void add_flux_divergence_inplace(double wght, torch::TensorList flux,
                                 torch::TensorList area,
                                 torch::Tensor const& vol, torch::Tensor& out) {
  // vol and area are 3D
  auto nx1 = vol.size(DIM1 - 1);
  auto nx2 = vol.size(DIM2 - 1);
  auto nx3 = vol.size(DIM3 - 1);

  // flux and out are 4D
  auto dflx =
      area[0].slice(DIM1 - 1, 1, nx1) * flux[0].slice(DIM1, 1, nx1) -
      area[0].slice(DIM1 - 1, 0, nx1 - 1) * flux[0].slice(DIM1, 0, nx1 - 1) +
      area[1].slice(DIM2 - 1, 1, nx2) * flux[1].slice(DIM2, 1, nx2) -
      area[1].slice(DIM2 - 1, 0, nx2 - 1) * flux[1].slice(DIM2, 0, nx2 - 1) +
      area[2].slice(DIM3 - 1, 1, nx3) * flux[2].slice(DIM3, 1, nx3) -
      area[2].slice(DIM3 - 1, 0, nx3 - 1) * flux[2].slice(DIM3, 0, nx3 - 1);

  out -= wght * dflx / vol;
}
