// torch
#include <torch/torch.h>

// exo3
#include <exo3/exo3.hpp>  // vec_raise_inplace

// snap
#include "eos.hpp"

enum {
  DIMC = 0,
  IDN = 0,
};

namespace canoe {
using namespace eos;

torch::Tensor eos_cons2prim_hydro_ideal(
    torch::Tensor &cons, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu, std::optional<torch::Tensor> rcv,
    int64_t ncloud, std::optional<torch::TensorList> cos_theta) {
  auto IVX = cons.size(DIMC) - 4;
  auto IVY = IVX + 1;
  auto IVZ = IVX + 2;
  auto IPR = IVX + 3;

  _check_dim1(IVX, ncloud, cons, gammad);

  _apply_conserved_limiter_inplace(cons);

  auto prim = torch::zeros_like(cons);

  prim[IDN] = cons.slice(DIMC, 0, IVX).sum(DIMC);
  prim.slice(DIMC, 1, IPR) = cons.slice(DIMC, 1, IPR) / prim[IDN];

  if (cos_theta.has_value()) {  // velocities are at an angle
    vec_raise_inplace(prim, cos_theta.value());
  }

  auto feps = torch::ones_like(cons[0]);

  if (rmu.has_value() && IVX > 1) {
    _check_dim2(IVX, rmu);
    auto primu =
        prim.slice(DIMC, 1, IVX - ncloud).unfold(DIMC, IVX - ncloud - 1, 1);
    feps += torch::matmul(primu, rmu.value().slice(DIMC, 0, IVX - ncloud - 1))
                .squeeze(DIMC) -
            prim.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);
  }

  auto fsig = torch::ones_like(cons[0]);

  if (rcv.has_value() && IVX > 1) {
    _check_dim2(IVX, rcv);
    auto primu = prim.slice(DIMC, 1, IVX).unfold(DIMC, IVX - 1, 1);
    fsig += torch::matmul(primu, rcv.value()).squeeze(DIMC);
  }

  auto ke = 0.5 * (prim[IVX] * cons[IVX] + prim[IVY] * cons[IVY] +
                   prim[IVZ] * cons[IVZ]);

  prim[IPR] = (gammad - 1) * (cons[IPR] - ke) * feps / fsig;

  _apply_primitive_limiter_inplace(prim);

  return prim;
}

torch::Tensor eos_prim2cons_hydro_ideal(
    torch::Tensor &prim, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu, std::optional<torch::Tensor> rcv,
    int64_t ncloud, std::optional<torch::TensorList> cos_theta) {
  auto IVX = prim.size(DIMC) - 4;
  auto IVY = IVX + 1;
  auto IVZ = IVX + 2;
  auto IPR = IVX + 3;

  _check_dim1(IVX, ncloud, prim, gammad);

  _apply_primitive_limiter_inplace(prim);

  auto cons = torch::zeros_like(prim);

  cons[IDN] = (1. - prim.slice(DIMC, 1, IVX).sum(DIMC)) * prim[IDN];
  cons.slice(DIMC, 1, IPR) = prim.slice(DIMC, 1, IPR) * prim[IDN];

  if (cos_theta.has_value()) {  // velocities are at an angle
    vec_lower_inplace(cons, cos_theta.value());
  }

  auto feps = torch::ones_like(cons[0]);

  if (rmu.has_value()) {
    _check_dim2(IVX, rmu);
    auto primu =
        prim.slice(DIMC, 1, IVX - ncloud).unfold(DIMC, IVX - ncloud - 1, 1);

    feps += torch::matmul(primu, rmu.value().slice(DIMC, 0, IVX - ncloud - 1))
                .squeeze(DIMC) -
            prim.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);
  }

  auto fsig = torch::ones_like(cons[0]);

  if (rcv.has_value()) {
    _check_dim2(IVX, rcv);
    auto primu = prim.slice(DIMC, 1, IVX).unfold(DIMC, IVX - 1, 1);
    fsig += torch::matmul(primu, rcv.value()).squeeze(DIMC);
  }

  auto ke = 0.5 * (prim[IVX] * cons[IVX] + prim[IVY] * cons[IVY] +
                   prim[IVZ] * cons[IVZ]);

  cons[IPR] = prim[IPR] * fsig / feps / (gammad - 1) + ke;

  _apply_conserved_limiter_inplace(cons);

  return cons;
}
}  // namespace canoe
