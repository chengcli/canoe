// torch
#include <torch/torch.h>

// snap
#include "hydro.hpp"

namespace canoe {

enum {
  DIM1 = 3,
  DIM2 = 2,
  DIM3 = 1,
  DIMC = 0,
};

torch::Tensor flux_divergence(torch::TensorList flux, torch::TensorList area,
                              torch::Tensor const& vol) {
  // vol and area are 3D
  auto nvar = flux[0].size(DIMC);
  auto nx1 = vol.size(DIM1);
  auto nx2 = vol.size(DIM2);
  auto nx3 = vol.size(DIM3);

  torch::Tensor dflx = torch::zeros({nvar, nx3, nx2, nx1}, vol.options());

  if (nx1 > 1) {
    dflx.slice(DIM1, 0, nx1 - 1) +=
        area[0].slice(DIM1, 1, nx1) * flux[0].slice(DIM1, 1, nx1) -
        area[0].slice(DIM1, 0, nx1 - 1) * flux[0].slice(DIM1, 0, nx1 - 1);
  }

  if (nx2 > 1) {
    dflx.slice(DIM2, 0, nx2 - 1) +=
        area[1].slice(DIM2, 1, nx2) * flux[1].slice(DIM2, 1, nx2) -
        area[1].slice(DIM2, 0, nx2 - 1) * flux[1].slice(DIM2, 0, nx2 - 1);
  }

  if (nx3 > 1) {
    dflx.slice(DIM3, 0, nx3 - 1) +=
        area[2].slice(DIM3, 1, nx3) * flux[2].slice(DIM3, 1, nx3) -
        area[2].slice(DIM3, 0, nx3 - 1) * flux[2].slice(DIM3, 0, nx3 - 1);
  }

  return dflx / vol;
}

}  // namespace canoe
