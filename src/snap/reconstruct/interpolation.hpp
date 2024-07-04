#pragma once

#define sqr(x) ((x) * (x))
#include <torch/torch.h>

namespace torch {

namespace cp5coeff {
Tensor c1m =
    tensor({-1. / 20., 9. / 20., 47. / 60., -13. / 60., 1. / 30.}, kFloat32);
Tensor c1p = flip(c1m, {0});
}  // namespace cp5coeff

inline Tensor interp_cp5m(Tensor const& phi, std::string pat) {
  return einsum(pat, {phi, cp5coeff::c1m});
};

inline Tensor interp_cp5p(Tensor const& phi, std::string pat) {
  return einsum(pat, {phi, cp5coeff::c1p});
};

namespace weno5coeff {
Tensor c1m = tensor({-1. / 6., 5. / 6., 1. / 3., 0., 0.}, kFloat32);
Tensor c1p = flip(c1m, {0});

Tensor c2m = tensor({0., 1. / 3., 5. / 6., -1. / 6., 0.}, kFloat32);
Tensor c2p = flip(c2m, {0});

Tensor c3m = tensor({0., 0., 11. / 6., -7. / 6., 1. / 3.}, kFloat32);
Tensor c3p = flip(c3m, {0});

Tensor c4m = tensor({1, -2, 1, 0, 0}, kFloat32);
Tensor c4p = flip(c4m, {0});

Tensor c5m = tensor({1, -4, 3, 0, 0}, kFloat32);
Tensor c5p = flip(c5m, {0});

Tensor c6m = tensor({0, 1, -2, 1, 0}, kFloat32);
Tensor c6p = flip(c6m, {0});

Tensor c7m = tensor({0, -1, 0, 1, 0}, kFloat32);
Tensor c7p = flip(c7m, {0});

Tensor c8m = tensor({0, 0, 1, -2, 1}, kFloat32);
Tensor c8p = flip(c8m, {0});

Tensor c9m = tensor({0, 0, 3, -4, 1}, kFloat32);
Tensor c9p = flip(c9m, {0});

}  // namespace weno5coeff

inline Tensor interp_weno5m(Tensor const& phi, std::string pat) {
  Tensor p1 = einsum(pat, {phi, weno5coeff::c1m});
  Tensor p2 = einsum(pat, {phi, weno5coeff::c2m});
  Tensor p3 = einsum(pat, {phi, weno5coeff::c3m});

  Tensor beta1 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c4m})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c5m}));
  Tensor beta2 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c6m})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c7m}));
  Tensor beta3 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c8m})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c9m}));

  Tensor alpha1 = 0.3 / sqr(beta1 + 1e-6);
  Tensor alpha2 = 0.6 / sqr(beta2 + 1e-6);
  Tensor alpha3 = 0.1 / sqr(beta3 + 1e-6);
  return (alpha1 * p1 + alpha2 * p2 + alpha3 * p3) / (alpha1 + alpha2 + alpha3);
}

inline Tensor interp_weno5p(Tensor const& phi, std::string pat) {
  Tensor p1 = einsum(pat, {phi, weno5coeff::c1p});
  Tensor p2 = einsum(pat, {phi, weno5coeff::c2p});
  Tensor p3 = einsum(pat, {phi, weno5coeff::c3p});

  Tensor beta1 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c4p})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c5p}));
  Tensor beta2 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c6p})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c7p}));
  Tensor beta3 = 13. / 12. * sqr(einsum(pat, {phi, weno5coeff::c8p})) +
                 1. / 4. * sqr(einsum(pat, {phi, weno5coeff::c9p}));

  Tensor alpha1 = 0.3 / sqr(beta1 + 1e-6);
  Tensor alpha2 = 0.6 / sqr(beta2 + 1e-6);
  Tensor alpha3 = 0.1 / sqr(beta3 + 1e-6);
  return (alpha1 * p1 + alpha2 * p2 + alpha3 * p3) / (alpha1 + alpha2 + alpha3);
}

}  // namespace torch
