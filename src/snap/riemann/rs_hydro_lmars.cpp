// torch
#include <torch/torch.h>

#define sqr(x) ((x) * (x))

enum {
  DIMC = 0,
  IDN = 0,
};

torch::Tensor rs_hydro_lmars(int8_t ivx, int8_t IVX, int8_t IPR,
                             const torch::Tensor &wl, const torch::Tensor &wr,
                             double gammad, std::optional<torch::Tensor> rmu,
                             std::optional<torch::Tensor> rcv) {
  int ivy = IVX + ((ivx - IVX) + 1) % 3;
  int ivz = IVX + ((ivx - IVX) + 2) % 3;
  int dir = ivx - IVX;

  auto fepsl = torch::ones_like(wl.slice(DIMC, 1, IVX));
  auto fepsr = torch::ones_like(wr.slice(DIMC, 1, IVX));

  auto fsigl = torch::ones_like(wl.slice(DIMC, 1, IVX));
  auto fsigr = torch::ones_like(wr.slice(DIMC, 1, IVX));

  if (rmu.has_value()) {
    fepsl += torch::einsum("nkji,n->kji",
                           {wl.slice(DIMC, 1, IVX), rmu.value() - 1.});
    fepsr += torch::einsum("nkji,n->kji",
                           {wr.slice(DIMC, 1, IVX), rmu.value() - 1.});
  }

  if (rcv.has_value()) {
    fsigl += torch::einsum("nkji,n->kji",
                           {wl.slice(DIMC, 1, IVX), rcv.value() - 1.});
    fsigr += torch::einsum("nkji,n->kji",
                           {wr.slice(DIMC, 1, IVX), rcv.value() - 1.});
  }

  auto kappal = 1.f / (gammad - 1.f) * fsigl / fepsl;
  auto kappar = 1.f / (gammad - 1.f) * fsigr / fepsr;

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

  auto fluxl = torch::zeros_like(wl);
  auto fluxr = torch::zeros_like(wr);

  auto rdryl = torch::ones_like(wl[IDN]);
  auto rdryr = torch::ones_like(wr[IDN]);

  rdryl -= wl.slice(DIMC, 1, IVX).sum(DIMC);
  rdryr -= wr.slice(DIMC, 1, IVX).sum(DIMC);

  fluxl[IDN] = ubar * wl[IDN] * rdryl;
  fluxl.slice(DIMC, 1, IVX) = ubar * wl[IDN] * wl.slice(DIMC, 1, IVX);

  fluxl[ivx] = ubar * wl[IDN] * wl[ivx] + pbar;
  fluxl[ivy] = ubar * wl[IDN] * wl[ivy];
  fluxl[ivz] = ubar * wl[IDN] * wl[ivz];
  fluxl[IPR] = ubar * wl[IDN] * hl;

  fluxr[IDN] = ubar * wr[IDN] * rdryr;
  fluxr.slice(DIMC, 1, IVX) = ubar * wr[IDN] * wr.slice(DIMC, 1, IVX);

  fluxr[ivx] = ubar * wr[IDN] * wr[ivx] + pbar;
  fluxr[ivy] = ubar * wr[IDN] * wr[ivy];
  fluxr[ivz] = ubar * wr[IDN] * wr[ivz];
  fluxr[IPR] = ubar * wr[IDN] * hr;

  auto ui = (ubar > 0).to(torch::kInt);

  return ui * fluxl + (1 - ui) * fluxr;
}
