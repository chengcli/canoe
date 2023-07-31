//! \file  roe_shallow_yz.cpp
//  \brief Roe's linearized Riemann solver for shallow water model

// C/C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <iostream>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          AthenaArray<Real> const &dxw) {
  int ivy = ivx == IVY ? IVZ : IVY;

  Real wli[NHYDRO], wri[NHYDRO], wave[3][3], speed[3];

  Real ubar, vbar, cbar, delh, delu, delv, hbar, a1, a2, a3;

  AthenaArray<Real> empty{};  // placeholder for unused electric/magnetic fields
#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of projection
  {
    switch (ivx) {
      case IVY:
        pmy_block->pcoord->PrimToLocal2(k, j, il, iu, empty, wl, wr, empty);
        break;
      case IVZ:
        pmy_block->pcoord->PrimToLocal3(k, j, il, iu, empty, wl, wr, empty);
        break;
    }
  }
#endif  // AFFINE or CUBED_SPHERE

  for (int i = il; i <= iu; ++i) {
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }

    ubar = (wli[ivx] * sqrt(wli[IDN]) + wri[ivx] * sqrt(wri[IDN])) /
           (sqrt(wli[IDN]) + sqrt(wri[IDN]));
    vbar = (wli[ivy] * sqrt(wli[IDN]) + wri[ivy] * sqrt(wri[IDN])) /
           (sqrt(wli[IDN]) + sqrt(wri[IDN]));
    cbar = sqrt(0.5 * (wli[IDN] + wri[IDN]));

    delh = wri[IDN] - wli[IDN];
    delu = wri[ivx] - wli[ivx];
    delv = wri[ivy] - wli[ivy];
    hbar = sqrt(wli[IDN] * wri[IDN]);

    a1 = 0.5 * (cbar * delh - hbar * delu) / cbar;
    a2 = hbar * delv;
    a3 = 0.5 * (cbar * delh + hbar * delu) / cbar;

    wave[0][0] = a1;
    wave[0][1] = a1 * (ubar - cbar);
    wave[0][2] = a1 * vbar;
    wave[1][0] = 0.;
    wave[1][1] = 0.;
    wave[1][2] = a2;
    wave[2][0] = a3;
    wave[2][1] = a3 * (ubar + cbar);
    wave[2][2] = a3 * vbar;

    speed[0] = fabs(ubar - cbar);
    speed[1] = fabs(ubar);
    speed[2] = fabs(ubar + cbar);

    flx(IDN, k, j, i) = 0.5 * (wli[IDN] * wli[ivx] + wri[IDN] * wri[ivx]);
    flx(ivx, k, j, i) =
        0.5 * (wli[IDN] * wli[ivx] * wli[ivx] + 0.5 * wli[IDN] * wli[IDN] +
               wri[IDN] * wri[ivx] * wri[ivx] + 0.5 * wri[IDN] * wri[IDN]);
    flx(ivy, k, j, i) =
        0.5 * (wli[IDN] * wli[ivx] * wli[ivy] + wri[IDN] * wri[ivx] * wri[ivy]);
    flx(IVX, k, j, i) = 0.;

    for (int r = 0; r < 3; ++r) {
      flx(IDN, k, j, i) -= 0.5 * speed[r] * wave[r][0];
      flx(ivx, k, j, i) -= 0.5 * speed[r] * wave[r][1];
      flx(ivy, k, j, i) -= 0.5 * speed[r] * wave[r][2];
    }
  }

#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of deprojection
  {
    // check if this is correct
    switch (ivx) {
      case IVY:
        pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, empty, empty, flx, empty,
                                         empty);
        break;
      case IVZ:
        pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, empty, empty, flx, empty,
                                         empty);
        break;
    }
  }
#endif  // AFFINE or CUBED_SPHERE

#if defined(CUBED_SPHERE)  // projection from contravariant fluxes to covariant
  {
    switch (ivx) {
      case IVY:
        for (int i = il; i <= iu; ++i) {
          Real x = tan(pmy_block->pcoord->x2f(j));
          Real y = tan(pmy_block->pcoord->x3v(k));
          Real C = sqrt(1. + x * x);
          Real D = sqrt(1. + y * y);
          Real cth = -x * y / C / D;

          // Extract local conserved quantities and fluxes
          const Real ty = flx(ivx, k, j, i);
          const Real tz = flx(ivy, k, j, i);

          // Transform fluxes
          flx(ivx, k, j, i) = ty + tz * cth;
          flx(ivy, k, j, i) = tz + ty * cth;
        }
        break;

      case IVZ:
        for (int i = il; i <= iu; ++i) {
          Real x = tan(pmy_block->pcoord->x2v(j));
          Real y = tan(pmy_block->pcoord->x3f(k));
          Real C = sqrt(1. + x * x);
          Real D = sqrt(1. + y * y);
          Real cth = -x * y / C / D;

          // Extract local conserved quantities and fluxes
          const Real ty = flx(ivx, k, j, i);
          const Real tz = flx(ivy, k, j, i);

          // Transform fluxes
          flx(ivx, k, j, i) = ty + tz * cth;
          flx(ivy, k, j, i) = tz + ty * cth;
        }
        break;
      default:
        break;
    }
  }

#endif
}
