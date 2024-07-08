// torch
#include <torch/torch.h>

enum {
  DIM1 = 2,
  IDN = 0,
};

namespace canoe {

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor> decom_obtain_anomaly(
    int64_t is, int64_t ie, torch::Tensor const& w, torch::Tensor const& dx1f,
    torch::Tensor const& grav) {
  auto sizes = w.sizes();
  auto nx1 = w.size(DIM1 + 1);
  auto wa = torch::zeros_like(w);

  auto IPR = w.size(0) - 1;

  auto psf =
      torch::zeros({w.size(0), w.size(1), w.size(2), nx1 + 1}, w.options());
  auto tsf = torch::zeros(psf.sizes(), w.options());

  // pressure anomaly
  auto RdTv = w[IPR].select(DIM1, ie) / w[IDN].select(DIM1, ie);
  auto rhogz =
      -grav * w[IDN].slice(DIM1, is, ie + 1) * dx1f.slice(DIM1, is, ie + 1);

  psf.select(DIM1, ie + 1) =
      w[IPR].select(DIM1, ie) * exp(-grav * dx1f[ie] / (2. * RdTv));
  psf.slice(DIM1, is + 1, ie + 1) = rhogz.slice(DIM1, is, ie);
  psf.select(DIM1, is) = psf.select(DIM1, ie + 1) + rhogz.sum(DIM1);
  psf.slice(DIM1, is, ie + 1) = psf.slice(DIM1, is, ie + 1).cumsum(DIM1);

  auto dpsf = psf.slice(DIM1, 0, nx1 - 1) - psf.slice(DIM1, 1, nx1);
  auto mask = dpsf.abs() < 1.e-5;

  auto psv = torch::where(
      mask, 0.5 * (psf.slice(DIM1, 0, nx1 - 1) + psf.slice(DIM1, 1, nx1)),
      dpsf / (psf.slice(DIM1, 0, nx1 - 1) / psf.slice(DIM1, 1, nx1)).log());

  //! virtual temperature anomaly
  auto tsv = psv / w[IDN];
  tsf.slice(DIM1, is + 1, ie) =
      0.5 * (tsv.slice(DIM1, is, ie - 1) + tsv.slice(DIM1, is + 1, ie));
  tsf.select(DIM1, ie + 1) = tsv.select(DIM1, ie);
  tsf.select(DIM1, is) = 2 * tsv.select(DIM1, is) - tsf.select(DIM1, is + 1);

  wa[IDN] = w[IPR] / w[IDN] - tsv;
  wa[IPR] -= psv;
  return {wa, psf, tsf};
}

void decom_apply_anomaly_inplace(torch::Tensor const& psf,
                                 torch::Tensor const& tsf, torch::Tensor& wl,
                                 torch::Tensor& wr) {
  auto IPR = wl.size(0) - 1;

  wl[IPR] += psf;
  torch::Tensor mask = wl[IPR] < 0;
  wl[IPR] = torch::where(mask, psf, wl[IPR]);

  wr[IPR] += psf;
  mask = wr[IPR] < 0;
  wr[IPR] = torch::where(mask, psf, wr[IPR]);

  wl[IDN] += tsf;
  mask = wl[IDN] < 0;
  wl[IDN] = wl[IPR] / torch::where(mask, tsf, wl[IDN]);

  wr[IDN] += tsf;
  mask = wr[IDN] < 0;
  wr[IDN] = wr[IPR] / torch::where(mask, tsf, wr[IDN]);
}

}  // namespace canoe
