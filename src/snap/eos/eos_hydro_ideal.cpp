// torch
#include <torch/torch.h>

// exo3
#include <exo3/exo3.hpp>  // vec_raise_inplace

// snap
#include "eos.hpp"

using namespace eos::internal;

enum {
  DIMC = 0,
  IDN = 0,
};

torch::Tensor eos_cons2prim_hydro_ideal(
    int64_t IVX, torch::Tensor &cons, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu, std::optional<torch::Tensor> rcv,
    int64_t ncloud, std::optional<torch::TensorList> cos_theta) {
  auto IVY = IVX + 1;
  auto IVZ = IVX + 2;
  auto IPR = IVX + 3;

  _apply_conserved_limiter_inplace(cons);

  auto prim = torch::zeros_like(cons);

  prim[IDN] = cons / cons.sum(DIMC);
  prim.slice(DIMC, IVX, IPR) = cons.slice(DIMC, IVX, IPR) / prim[IDN];

  auto feps = torch::ones_like(cons[0]);

  if (rmu.has_value()) {
    feps += torch::einsum("nkji,n->kji",
                          {prim.slice(DIMC, 1, IVX - ncloud),
                           rmu.value().slice(DIMC, 1, IVX - ncloud)}) -
            prim.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);
  }

  auto fsig = torch::ones_like(cons[0]);

  if (rcv.has_value()) {
    fsig +=
        torch::einsum("nkji,n->kji", {prim.slice(DIMC, 1, IVX), rcv.value()});
  }

  if (cos_theta.has_value()) {  // velocities are at an angle
    vec_raise_inplace(prim, cos_theta.value());
  }

  auto ke = 0.5 * (prim[IVX] * cons[IVX] + prim[IVY] * cons[IVY] +
                   prim[IVZ] * cons[IVZ]);

  prim[IPR] = (gammad - 1) * (cons[IPR] - ke) * feps / fsig;

  _apply_primitive_limiter_inplace(prim);

  return prim;
}

torch::Tensor eos_prim2cons_hydro_ideal(
    int64_t IVX, torch::Tensor &prim, const torch::Tensor &gammad,
    std::optional<torch::Tensor> rmu, std::optional<torch::Tensor> rcv,
    int64_t ncloud, std::optional<torch::TensorList> cos_theta) {
  auto IVY = IVX + 1;
  auto IVZ = IVX + 2;
  auto IPR = IVX + 3;

  _apply_primitive_limiter_inplace(prim);

  auto cons = torch::zeros_like(prim);

  cons[IDN] = 1. - prim.slice(DIMC, 1, IVX).sum(DIMC);
  cons.slice(DIMC, 1, IVX) = prim.slice(DIMC, 1, IVX) * prim[IDN];
  cons.slice(DIMC, IVX, IPR) = prim.slice(DIMC, IVX, IPR) * prim[IDN];

  auto feps = torch::ones_like(cons[0]);

  if (rmu.has_value()) {
    feps += torch::einsum("nkji,n->kji",
                          {prim.slice(DIMC, 1, IVX - ncloud),
                           rmu.value().slice(DIMC, 1, IVX - ncloud)}) -
            prim.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);
  }

  auto fsig = torch::ones_like(cons[0]);

  if (rcv.has_value()) {
    fsig +=
        torch::einsum("nkji,n->kji", {prim.slice(DIMC, 1, IVX), rcv.value()});
  }

  if (cos_theta.has_value()) {  // velocities are at an angle
    vec_lower_inplace(cons, cos_theta.value());
  }

  auto ke = 0.5 * (prim[IVX] * cons[IVX] + prim[IVY] * cons[IVY] +
                   prim[IVZ] * cons[IVZ]);

  cons[IPR] = prim[IPR] * fsig / feps / (gammad - 1) + ke;

  _apply_conserved_limiter_inplace(cons);

  return cons;
}
