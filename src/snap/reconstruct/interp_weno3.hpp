#ifndef INTERP_WENO3_HPP_
#define INTERP_WENO3_HPP_
#include "../defs.hpp"

namespace Weno3Coeff {
double a1[2] = {RAT1 * RAT1 / (1 + RAT1), 1 / (1 + RAT1)};

double a2[2] = {(2 + RAT1) / (1 + RAT1), -1 / (RAT1 * (1 + RAT1))};

double wght[2] = {RAT1 * (1 + RAT1) / (1 + RAT1 + RAT1 * RAT1),
                  1 / (1 + RAT1 + RAT1 * RAT1)};

double beta1[2][3] = {{11. / 4., -3. / 2., -1. / 4.}, {1. / 2., 1. / 2., 0.}};

double beta2[2][3] = {{3. / 4., -11. / 2., 23. / 4.},
                      {-1. / 2., 13. / 2., -9.}};

double c1[3] = {
    RAT1 * RAT1 * RAT1 / (1 + RAT1 + RAT1 * RAT1),
    (2 + 2 * RAT1 + RAT1 * RAT1) / ((1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
    -1 / (RAT1 * (1 + RAT1) * (1 + RAT1 + RAT1 * RAT1))};
}  // namespace Weno3Coeff

inline Real interp_weno3x1(Real phim1, Real phi, Real phip1) {
  using namespace Weno3Coeff;
  Real p1 = a1[0] * phim1 + a1[1] * phi;
  Real p2 = a2[0] * phi + a2[1] * phip1;

  Real eps = RAT1 - 1.;
  Real pp1[3] = {phim1 * phim1, phim1 * phi, phi * phi};
  Real b1 =
      (phim1 - phi) * (phim1 - phi) + (3 * phim1 - phi) * (phim1 - phi) * eps;
  for (int i = 0; i < 3; ++i)
    b1 += (beta1[0][i] + beta1[1][i] * eps) * eps * eps * pp1[i];

  Real pp2[3] = {phi * phi, phi * phip1, phip1 * phip1};
  Real b2 =
      (phi - phip1) * (phi - phip1) - (phi - phip1) * (phi - 3 * phip1) * eps;
  for (int i = 0; i < 3; ++i)
    b2 += (beta2[0][i] + beta2[1][i] * eps) * eps * eps * pp2[i];

  Real w1 = wght[0] / ((b1 + 1.E-10) * (b1 + 1.E-10));
  Real w2 = wght[1] / ((b2 + 1.E-10) * (b2 + 1.E-10));

  return (w1 * p1 + w2 * p2) / (w1 + w2);
}

inline Real interp_cp3x1(Real phim1, Real phi, Real phip1) {
  using namespace Weno3Coeff;
  return c1[0] * phim1 + c1[1] * phi + c1[2] * phip1;
}

#endif
