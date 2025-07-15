//! \file lmars.cpp
//  \brief LMARS based on Chen's 2013 paper

// C/C++
#include <cstring>

// Athena++
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include <impl.hpp>
#include <interface/eos.hpp>

// exo3
#include <exo3/exo3.hpp>

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx - IVX) + 1) % 3;
  int ivz = IVX + ((ivx - IVX) + 2) % 3;
  int dir = ivx - IVX;

  MeshBlock *pmb = pmy_block;

  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pmb->peos->GetGamma();
  Real wli[NHYDRO], wri[NHYDRO];

  auto peos = pmb->pimpl->peos;
  for (int i = il; i <= iu; ++i) {
    // copy local variables
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }

    // correction for gamma
    Real kappal = 1. / (get_gamma(peos, wli) - 1.);
    Real kappar = 1. / (get_gamma(peos, wri) - 1.);

    // enthalpy
    // FIXME: m should be at cell interface
    auto vl = vec_lower<Real>(wli, pmb->pcoord->m.at(k, j, i));
    hl = wli[IPR] / wli[IDN] * (kappal + 1.) +
         0.5 * (wli[ivx] * vl[ivx] + wli[ivy] * vl[ivy] + wli[ivz] * vl[ivz]);

    auto vr = vec_lower<Real>(wri, pmb->pcoord->m.at(k, j, i));
    hr = wri[IPR] / wri[IDN] * (kappar + 1.) +
         0.5 * (wri[ivx] * vr[ivx] + wri[ivy] * vr[ivy] + wri[ivz] * vr[ivz]);

    rhobar = 0.5 * (wli[IDN] + wri[IDN]);
    cbar = sqrt(0.5 * (1. + (1. / kappar + 1. / kappal) / 2.) *
                (wli[IPR] + wri[IPR]) / rhobar);
    pbar = 0.5 * (wli[IPR] + wri[IPR]) +
           0.5 * (rhobar * cbar) * (wli[ivx] - wri[ivx]);
    ubar = 0.5 * (wli[ivx] + wri[ivx]) +
           0.5 / (rhobar * cbar) * (wli[IPR] - wri[IPR]);

    Real rdl = 1., rdr = 1.;
    for (int n = 1; n < IVX; ++n) {
      rdl -= wli[n];
      rdr -= wri[n];
    }

    if (ubar > 0.) {
      flx(IDN, k, j, i) = ubar * wli[IDN] * rdl;
      for (int n = 1; n < IVX; ++n) flx(n, k, j, i) = ubar * wli[IDN] * wli[n];
      flx(ivx, k, j, i) = ubar * wli[IDN] * wli[ivx] + pbar;
      flx(ivy, k, j, i) = ubar * wli[IDN] * wli[ivy];
      flx(ivz, k, j, i) = ubar * wli[IDN] * wli[ivz];
      flx(IEN, k, j, i) = ubar * wli[IDN] * hl;
    } else {
      flx(IDN, k, j, i) = ubar * wri[IDN] * rdr;
      for (int n = 1; n < IVX; ++n) flx(n, k, j, i) = ubar * wri[IDN] * wri[n];
      flx(ivx, k, j, i) = ubar * wri[IDN] * wri[ivx] + pbar;
      flx(ivy, k, j, i) = ubar * wri[IDN] * wri[ivy];
      flx(ivz, k, j, i) = ubar * wri[IDN] * wri[ivz];
      flx(IEN, k, j, i) = ubar * wri[IDN] * hr;
    }
  }
}
