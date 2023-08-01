//! \file  lmars_shallow_yz.cpp
//  \brief

// C/C++
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <iostream>

// Athena++
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
  auto pcoord = pmy_block->pcoord;

  Real wli[NHYDRO], wri[NHYDRO];

  Real ubar, cbar, hbar;

  AthenaArray<Real> empty{};  // placeholder for unused electric/magnetic fields
#if defined(AFFINE) or defined(CUBED_SPHERE)  // need of projection
  if (ivx == IVY) {
    pcoord->PrimToLocal2(k, j, il, iu, empty, wl, wr, empty);
  } else if (ivx == IVZ) {
    pcoord->PrimToLocal3(k, j, il, iu, empty, wl, wr, empty);
  }
#endif  // AFFINE or CUBED_SPHERE

  for (int i = il; i <= iu; ++i) {
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }

    cbar = sqrt(0.5 * (wli[IDN] + wri[IDN]));
    hbar = 0.5 * (wli[IDN] + wri[IDN]) + 0.5 * cbar * (wli[ivx] - wri[ivx]);
    ubar = 0.5 * (wli[ivx] + wri[ivx]) + 0.5 / cbar * (wli[IDN] - wri[IDN]);

    if (ubar > 0.) {
      flx(IDN, k, j, i) = ubar * wli[IDN];
      flx(ivx, k, j, i) = ubar * wli[IDN] * wli[ivx] + 0.5 * hbar * hbar;
      flx(ivy, k, j, i) = ubar * wli[IDN] * wli[ivy];
    } else {
      flx(IDN, k, j, i) = ubar * wri[IDN];
      flx(ivx, k, j, i) = ubar * wri[IDN] * wri[ivx] + 0.5 * hbar * hbar;
      flx(ivy, k, j, i) = ubar * wri[IDN] * wri[ivy];
    }
    flx(IVX, k, j, i) = 0.;
  }

#if defined(AFFINE) or defined(CUBED_SPHERE)  // need of deprojection
  if (ivx == IVY) {
    pcoord->FluxToGlobal2(k, j, il, iu, empty, empty, flx, empty, empty);
  } else if (ivx == IVZ) {
    pcoord->FluxToGlobal3(k, j, il, iu, empty, empty, flx, empty, empty);
  }
#endif  // AFFINE or CUBED_SPHERE

#if defined(CUBED_SPHERE)  // projection from contravariant fluxes to covariant
  Real x, y;

  if (ivx == IVY) {
    x = tan(pcoord->x2f(j));
    y = tan(pcoord->x3v(k));
  } else if (ivx == IVZ) {  // IVZ
    x = tan(pcoord->x2v(j));
    y = tan(pcoord->x3f(k));
  }

  Real C = sqrt(1. + x * x);
  Real D = sqrt(1. + y * y);
  Real cth = -x * y / C / D;

  for (int i = il; i <= iu; ++i) {
    // Extract local conserved quantities and fluxes
    Real ty = flx(ivx, k, j, i);
    Real tz = flx(ivy, k, j, i);

    // Transform fluxes
    flx(ivx, k, j, i) = ty + tz * cth;
    flx(ivy, k, j, i) = tz + ty * cth;
  }
#endif  // CUBED_SPHERE
}
