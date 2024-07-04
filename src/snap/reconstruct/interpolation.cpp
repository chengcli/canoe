#include "interpolation.hpp"

#define sqr(x) ((x) * (x))

using namespace torch;

namespace cp5coeff {
Tensor c1m = tensor({-1. / 20., 9. / 20., 47. / 60., -13. / 60., 1. / 30.},
                    dtype(kFloat32));
}  // namespace cp5coeff

namespace weno5coeff {
Tensor c1m = tensor({-1. / 6., 5. / 6., 1. / 3., 0., 0.}, dtype(kFloat32));
Tensor c2m = tensor({0., 1. / 3., 5. / 6., -1. / 6., 0.}, dtype(kFloat32));
Tensor c3m = tensor({0., 0., 11. / 6., -7. / 6., 1. / 3.}, dtype(kFloat32));
Tensor c4m = tensor({1, -2, 1, 0, 0}, dtype(kFloat32));
Tensor c5m = tensor({1, -4, 3, 0, 0}, dtype(kFloat32));
Tensor c6m = tensor({0, 1, -2, 1, 0}, dtype(kFloat32));
Tensor c7m = tensor({0, -1, 0, 1, 0}, dtype(kFloat32));
Tensor c8m = tensor({0, 0, 1, -2, 1}, dtype(kFloat32));
Tensor c9m = tensor({0, 0, 3, -4, 1}, dtype(kFloat32));
}  // namespace weno5coeff

Center5Interp::Center5Interp(c10::DeviceType dtype) {
  cm_ = {cp5coeff::c1m};
  cp_ = cm_;
  for (auto& c : cp_) c = flip(c, {0});

  ToDevice(dtype);
}

Tensor Center5Interp::left(Tensor const& phi, std::string pat) {
  return einsum(pat, {phi, cm_[0]});
}

Tensor Center5Interp::right(Tensor const& phi, std::string pat) {
  return einsum(pat, {phi, cp_[0]});
}

Weno5Interp::Weno5Interp(c10::DeviceType dtype) {
  cm_ = {weno5coeff::c1m, weno5coeff::c2m, weno5coeff::c3m,
         weno5coeff::c4m, weno5coeff::c5m, weno5coeff::c6m,
         weno5coeff::c7m, weno5coeff::c8m, weno5coeff::c9m};
  cp_ = cm_;
  for (auto& c : cp_) c = flip(c, {0});

  ToDevice(dtype);
}

Tensor Weno5Interp::left(Tensor const& phi, std::string pat) {
  Tensor p1 = einsum(pat, {phi, cm_[0]});
  Tensor p2 = einsum(pat, {phi, cm_[1]});
  Tensor p3 = einsum(pat, {phi, cm_[2]});

  Tensor beta1 = 13. / 12. * sqr(einsum(pat, {phi, cm_[3]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cm_[4]}));
  Tensor beta2 = 13. / 12. * sqr(einsum(pat, {phi, cm_[5]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cm_[6]}));
  Tensor beta3 = 13. / 12. * sqr(einsum(pat, {phi, cm_[7]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cm_[8]}));

  Tensor alpha1 = 0.3 / sqr(beta1 + 1e-6);
  Tensor alpha2 = 0.6 / sqr(beta2 + 1e-6);
  Tensor alpha3 = 0.1 / sqr(beta3 + 1e-6);
  return (alpha1 * p1 + alpha2 * p2 + alpha3 * p3) / (alpha1 + alpha2 + alpha3);
}

Tensor Weno5Interp::right(Tensor const& phi, std::string pat) {
  Tensor p1 = einsum(pat, {phi, cp_[0]});
  Tensor p2 = einsum(pat, {phi, cp_[1]});
  Tensor p3 = einsum(pat, {phi, cp_[2]});

  Tensor beta1 = 13. / 12. * sqr(einsum(pat, {phi, cp_[3]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cp_[4]}));
  Tensor beta2 = 13. / 12. * sqr(einsum(pat, {phi, cp_[5]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cp_[6]}));
  Tensor beta3 = 13. / 12. * sqr(einsum(pat, {phi, cp_[7]})) +
                 1. / 4. * sqr(einsum(pat, {phi, cp_[8]}));

  Tensor alpha1 = 0.3 / sqr(beta1 + 1e-6);
  Tensor alpha2 = 0.6 / sqr(beta2 + 1e-6);
  Tensor alpha3 = 0.1 / sqr(beta3 + 1e-6);
  return (alpha1 * p1 + alpha2 * p2 + alpha3 * p3) / (alpha1 + alpha2 + alpha3);
}
