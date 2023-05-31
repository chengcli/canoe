#ifndef SRC_SNAP_RECONSTRUCT_INTERP_WENO5_HPP_
#define SRC_SNAP_RECONSTRUCT_INTERP_WENO5_HPP_

// athena
#include <athena/defs.hpp>

namespace Weno5Coeff {
double a1[3] = {
    -RAT1 * RAT1 * RAT1 * RAT1 * RAT1 / ((1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
    RAT1* RAT1*(1 + 2 * RAT1 + 2 * RAT1 * RAT1) /
        ((1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
    1 / (1 + RAT1 + RAT1 * RAT1)};

double a2[3] = {
    RAT1 * RAT1 * RAT1 / (1 + RAT1 + RAT1 * RAT1),
    (2 + 2 * RAT1 + RAT1 * RAT1) / ((1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
    -1 / (RAT1 * (1 + RAT1) * (1 + RAT1 + RAT1 * RAT1))};

double a3[3] = {(3 + 4 * RAT1 + 3 * RAT1 * RAT1 + RAT1 * RAT1 * RAT1) /
                    ((1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
                -(1 + 3 * RAT1 + 2 * RAT1 * RAT1 + RAT1 * RAT1 * RAT1) /
                    (RAT1 * RAT1 * (1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)),
                1 / (RAT1 * RAT1 * RAT1 * (1 + RAT1 + RAT1 * RAT1))};

double wght[3] = {
    (RAT1 * RAT1 * RAT1 * RAT1 * (1 + RAT1 + RAT1 * RAT1)) /
        ((1 + RAT1 * RAT1) * (1 + RAT1 + RAT1 * RAT1 + RAT1 * RAT1 * RAT1 +
                              RAT1 * RAT1 * RAT1 * RAT1)),
    (RAT1 * (1 + RAT1) * (1 + RAT1 + RAT1 * RAT1)) /
        ((1 + RAT1 * RAT1) * (1 + RAT1 + RAT1 * RAT1 + RAT1 * RAT1 * RAT1 +
                              RAT1 * RAT1 * RAT1 * RAT1)),
    1 / ((1 + RAT1 * RAT1) * (1 + RAT1 + RAT1 * RAT1 + RAT1 * RAT1 * RAT1 +
                              RAT1 * RAT1 * RAT1 * RAT1))};

double beta1[3][6] = {
    {4. / 3., -19. / 3., 11. / 3., 25. / 3., -31. / 3., 10. / 3.},
    {35. / 3., -132. / 3., 59. / 3., 124. / 3., -104. / 3., 6.},
    {1592. / 36., -4631. / 36., 1483. / 36., 2993. / 36., -1463. / 36.,
     62. / 36.}};

double beta2[3][6] = {
    {4. / 3., -13. / 3., 5. / 3., 13. / 3., -13. / 3., 4. / 3.},
    {13. / 3., -23. / 3., 0., 0., 23. / 3., -13. / 3.},
    {140. / 36., -5. / 36., -59. / 36., -55. / 36., -281. / 36., 296. / 36.}};

double beta3[3][6] = {
    {10. / 3., -31. / 3., 11. / 3., 25. / 3., -19. / 3., 4. / 3.},
    {-6., 104. / 3., -59. / 3., -124. / 3., 132. / 3., -35. / 3.},
    {278. / 36., -2711. / 36., 2191. / 36., 4481. / 36., -6215. / 36.,
     2012. / 36.}};

double c1[5] = {-3. / 60., 27. / 60., 47. / 60., -13. / 60., 2. / 60.};
double c2[5] = {-33. / 120., 147. / 120., -133. / 120., 107. / 120.,
                -28. / 120.};
double c3[5] = {-387. / 720., 513. / 720., 673. / 720., -1487. / 720.,
                628. / 720.};
double c4[5] = {-513. / 1440., -513. / 1440., -673. / 1440., 5087. / 1440.,
                -3328. / 1440.};
}  // namespace Weno5Coeff

#ifdef STRETCHED_GRID

inline Real interp_weno5(Real phim2, Real phim1, Real phi, Real phip1,
                         Real phip2) {
  using Weno5Coeff::*;
  Real p1 = 0, p2 = 0, p3 = 0;
  Real phis[5] = {phim2, phim1, phi, phip1, phip2};
  for (int i = 0; i < 3; ++i) {
    p1 += a1[i] * phis[i];
    p2 += a2[i] * phis[1 + i];
    p3 += a3[i] * phis[2 + i];
  }

  Real b1 = 0, b2 = 0, b3 = 0;
  Real eps = RAT1 - 1.;
  Real pp1[6] = {phim2 * phim2, phim2 * phim1, phim2 * phi,
                 phim1 * phim1, phim1 * phi,   phi * phi};
  Real pp2[6] = {phim1 * phim1, phim1 * phi, phim1 * phip1,
                 phi * phi,     phi * phip1, phip1 * phip1};
  Real pp3[6] = {phi * phi,     phi * phip1,   phi * phip2,
                 phip1 * phip1, phip1 * phip2, phip2 * phip2};

  for (int i = 0; i < 6; ++i) {
    b1 += (beta1[0][i] + beta1[1][i] * eps + beta1[2][i] * eps * eps) * pp1[i];
    b2 += (beta2[0][i] + beta2[1][i] * eps + beta2[2][i] * eps * eps) * pp2[i];
    b3 += (beta3[0][i] + beta3[1][i] * eps + beta3[2][i] * eps * eps) * pp3[i];
  }

  Real a1 = wght[0] / ((b1 + 1.E-10) * (b1 + 1.E-10));
  Real a2 = wght[1] / ((b2 + 1.E-10) * (b2 + 1.E-10));
  Real a3 = wght[2] / ((b3 + 1.E-10) * (b3 + 1.E-10));

  return (a1 * p1 + a2 * p2 + a3 * p3) / (a1 + a2 + a3);
}

inline Real interp_cp5(Real phim2, Real phim1, Real phi, Real phip1,
                       Real phip2) {
  using Weno5Coeff::*;
  Real phis[5] = {phim2, phim1, phi, phip1, phip2};
  Real result = 0.;
  Real eps = RAT1 - 1.;
  for (int i = 0; i < 5; ++i)
    result +=
        (c1[i] + c2[i] * eps + c3[i] * eps * eps + c4[i] * eps * eps * eps) *
        phis[i];
  return result;
}

#endif  // STRETCHED_GRID

#endif  // SRC_SNAP_RECONSTRUCT_INTERP_WENO5_HPP_
