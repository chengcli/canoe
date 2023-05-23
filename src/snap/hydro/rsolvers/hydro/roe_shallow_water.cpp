//! \file  roe_shallow_water.cpp
//  \brief Roe's linearized Riemann solver for shallow water model

// C/C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <iostream>

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../hydro.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          AthenaArray<Real> const &dxw) {
  int ivy = ivx == IVX ? IVY : IVX;

  Real wli[NHYDRO], wri[NHYDRO], wave[3][3], speed[3];

  Real ubar, vbar, cbar, delh, delu, delv, hbar, a1, a2, a3;

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
    flx(IVZ, k, j, i) = 0.;

    for (int r = 0; r < 3; ++r) {
      flx(IDN, k, j, i) -= 0.5 * speed[r] * wave[r][0];
      flx(ivx, k, j, i) -= 0.5 * speed[r] * wave[r][1];
      flx(ivy, k, j, i) -= 0.5 * speed[r] * wave[r][2];
    }
  }
}
