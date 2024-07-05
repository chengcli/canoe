// torch
#include <torch/torch.h>

// snap
#include "riemann.hpp"

#define sqr(x) ((x) * (x))

enum {
  DIMC = 0,
  IDN = 0,
};

using namespace riemann::internal;

torch::Tensor rs_hydro_lmars(int64_t dim, const torch::Tensor &wl,
                             const torch::Tensor &wr,
                             const torch::Tensor &gammad,
                             std::optional<torch::Tensor> rmu,
                             std::optional<torch::Tensor> rcv, int64_t ncloud,
                             std::optional<torch::TensorList> cos_theta) {
  auto IVX = wl.size(DIMC) - 4;
  auto IPR = IVX + 3;  // pressure index

  // dim, ivx
  // 3, IVX
  // 2, IVX + 1
  // 1, IVX + 2
  auto ivx = IPR - dim;
  auto ivy = IVX + ((ivx - IVX) + 1) % 3;
  auto ivz = IVX + ((ivx - IVX) + 2) % 3;

  auto fepsl = torch::ones_like(wl[0]);
  auto fepsr = torch::ones_like(wr[0]);

  if (rmu.has_value()) {
    auto wlu =
        wl.slice(DIMC, 1, IVX - ncloud).unfold(DIMC, IVX - ncloud - 1, 1);
    fepsl += torch::matmul(wlu, rmu.value().slice(DIMC, 0, IVX - ncloud - 1))
                 .squeeze(DIMC) -
             wl.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);

    auto wru =
        wr.slice(DIMC, 1, IVX - ncloud).unfold(DIMC, IVX - ncloud - 1, 1);
    fepsr += torch::matmul(wru, rmu.value().slice(DIMC, 0, IVX - ncloud - 1))
                 .squeeze(DIMC) -
             wr.slice(DIMC, IVX - ncloud, IVX).sum(DIMC);
  }

  auto fsigl = torch::ones_like(wl[0]);
  auto fsigr = torch::ones_like(wr[0]);

  if (rcv.has_value()) {
    auto wlu = wl.slice(DIMC, 1, IVX).unfold(DIMC, IVX - 1, 1);
    fsigl += torch::matmul(wlu, rcv.value()).squeeze(DIMC);

    auto wru = wr.slice(DIMC, 1, IVX).unfold(DIMC, IVX - 1, 1);
    fsigr += torch::matmul(wru, rcv.value()).squeeze(DIMC);
  }

  auto kappal = 1.f / (gammad - 1.f) * fsigl / fepsl;
  auto kappar = 1.f / (gammad - 1.f) * fsigr / fepsr;

  if (cos_theta.has_value()) {  // velocities are at an angle
    _prim2local_inplace(ivx, wl, cos_theta.value());
  }

  auto kel = 0.5 * (sqr(wl[ivx]) + sqr(wl[ivy]) + sqr(wl[ivz]));
  auto ker = 0.5 * (sqr(wr[ivx]) + sqr(wr[ivy]) + sqr(wr[ivz]));

  // enthalpy
  auto hl = wl[IPR] / wl[IDN] * (kappal + 1.) + kel;
  auto hr = wr[IPR] / wr[IDN] * (kappar + 1.) + ker;

  auto rhobar = 0.5 * (wl[IDN] + wr[IDN]);
  auto cbar = sqrt(0.5 * (1. + (1. / kappar + 1. / kappal) / 2.) *
                   (wl[IPR] + wr[IPR]) / rhobar);
  auto pbar =
      0.5 * (wl[IPR] + wr[IPR]) + 0.5 * (rhobar * cbar) * (wl[ivx] - wr[ivx]);
  auto ubar =
      0.5 * (wl[ivx] + wr[ivx]) + 0.5 / (rhobar * cbar) * (wl[IPR] - wr[IPR]);

  // left/right flux
  auto fluxl = torch::zeros_like(wl);
  auto fluxr = torch::zeros_like(wr);

  fluxl[IDN] = ubar * wl[IDN] *
               (torch::ones_like(wl[IDN]) - wl.slice(DIMC, 1, IVX).sum(DIMC));
  fluxl.slice(DIMC, 1, IVX) = ubar * wl[IDN] * wl.slice(DIMC, 1, IVX);

  fluxl[ivx] = ubar * wl[IDN] * wl[ivx] + pbar;
  fluxl[ivy] = ubar * wl[IDN] * wl[ivy];
  fluxl[ivz] = ubar * wl[IDN] * wl[ivz];
  fluxl[IPR] = ubar * wl[IDN] * hl;

  fluxr[IDN] = ubar * wr[IDN] *
               (torch::ones_like(wr[IDN]) - wr.slice(DIMC, 1, IVX).sum(DIMC));
  fluxr.slice(DIMC, 1, IVX) = ubar * wr[IDN] * wr.slice(DIMC, 1, IVX);

  fluxr[ivx] = ubar * wr[IDN] * wr[ivx] + pbar;
  fluxr[ivy] = ubar * wr[IDN] * wr[ivy];
  fluxr[ivz] = ubar * wr[IDN] * wr[ivz];
  fluxr[IPR] = ubar * wr[IDN] * hr;

  auto ui = (ubar > 0).to(torch::kInt);
  auto flx = ui * fluxl + (1 - ui) * fluxr;

  if (cos_theta.has_value()) {  // velocities are at an angle
    _flux2global_inplace(ivx, flx, cos_theta.value());
  }

  return flx;
}
