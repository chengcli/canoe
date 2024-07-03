#pragma once

#define sqr(x) ((x) * (x))
#include <torch/torch.h>

namespace torch {

namespace weno5coeff {
Tensor c1 = tensor({-1. / 6., 5. / 6., 1. / 3., 0., 0.}, kFloat32);
Tensor c2 = tensor({0., 1. / 3., 5. / 6., -1. / 6., 0.}, kFloat32);
Tensor c3 = tensor({0., 0., 11. / 6., -7. / 6., 1. / 3.}, kFloat32);

Tensor c4 = tensor({1, -2, 1, 0, 0}, kFloat32);
Tensor c5 = tensor({1, -4, 3, 0, 0}, kFloat32);

Tensor c6 = tensor({0, 1, -2, 1, 0}, kFloat32);
Tensor c7 = tensor({0, -1, 0, 1, 0}, kFloat32);

Tensor c8 = tensor({0, 0, 1, -2, 1}, kFloat32);
Tensor c9 = tensor({0, 0, 3, -4, 1}, kFloat32);

Tensor c10 = tensor({0.3, 0.6, 0.1}, kFloat32);
}  // namespace weno5coeff

inline float interp_weno5(Tensor const& phi) {
  Tensor p = zeros(3, kFloat32);
  p[0] = dot(weno5coeff::c1, phi);
  p[1] = dot(weno5coeff::c2, phi);
  p[2] = dot(weno5coeff::c3, phi);

  Tensor beta = zeros(3, kFloat32);
  beta[0] = 13. / 12. * sqr(dot(weno5coeff::c4, phi)) +
            1. / 4. * sqr(dot(weno5coeff::c5, phi));
  beta[1] = 13. / 12. * sqr(dot(weno5coeff::c6, phi)) +
            1. / 4. * sqr(dot(weno5coeff::c7, phi));
  beta[2] = 13. / 12. * sqr(dot(weno5coeff::c8, phi)) +
            1. / 4. * sqr(dot(weno5coeff::c9, phi));

  Tensor alpha = weno5coeff::c10 / sqr(beta + 1e-10);
  return (dot(alpha, p) / sum(alpha)).item<float>();
}

}  // namespace torch
