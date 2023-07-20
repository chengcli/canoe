//! \file  lmars_shallow_yz.cpp
//  \brief

// C/C++
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <iostream>

// Athena++
#include <athena.hpp>
#include <coordinates/coordinates.hpp>
#include <hydro/hydro.hpp>
#include <mesh/mesh.hpp>

// canoe
#include <configure.hpp>

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          AthenaArray<Real> const &dxw) {
  int ivy = ivx == IVY ? IVZ : IVY;

  Real wli[NHYDRO], wri[NHYDRO];

  Real ubar, cbar, hbar;

  AthenaArray<Real> empty{};  // placeholder for unused electric/magnetic fields
#ifdef AFFINE                 // need of projection
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
#endif  // AFFINE

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
      flx(ivx, k, j, i) = ubar * wri[IDN] * wri[ivx] * 0.5 * hbar * hbar;
      flx(ivy, k, j, i) = ubar * wri[IDN] * wri[ivy];
    }
    flx(IVX, k, j, i) = 0.;
  }

#ifdef AFFINE  // need of deprojection
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
#endif  // AFFINE
}
