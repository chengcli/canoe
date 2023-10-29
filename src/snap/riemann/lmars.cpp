//! \file lmars.cpp
//  \brief LMARS based on Chen's 2013 paper

// C/C++
#include <cstring>

// Athena++
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// climath
#include <climath/core.h>  // sqr

// canoe
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// microphysics
#include <microphysics/microphysics.hpp>

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          const AthenaArray<Real> &dxw) {
  int ivy = IVX + ((ivx - IVX) + 1) % 3;
  int ivz = IVX + ((ivx - IVX) + 2) % 3;
  int dir = ivx - IVX;

  auto pthermo = Thermodynamics::GetInstance();
  auto pmicro = pmy_block->pimpl->pmicro;

  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pmy_block->peos->GetGamma();
  Real wli[NHYDRO], wri[NHYDRO];

  for (int i = il; i <= iu; ++i) {
    // copy local variables
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }
    // correction for gamma
    // left
    Real fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += wli[n] * (pthermo->GetCvRatioMass(n) - 1.);
      feps += wli[n] * (1. / pthermo->GetMuRatio(n) - 1.);
    }
    Real kappal = 1. / (gamma - 1.) * fsig / feps;

    // right
    fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += wri[n] * (pthermo->GetCvRatioMass(n) - 1.);
      feps += wri[n] * (1. / pthermo->GetMuRatio(n) - 1.);
    }
    Real kappar = 1. / (gamma - 1.) * fsig / feps;

    // enthalpy
    hl = wli[IPR] / wli[IDN] * (kappal + 1.) +
         0.5 * (sqr(wli[ivx]) + sqr(wli[ivy]) + sqr(wli[ivz]));
    hr = wri[IPR] / wri[IDN] * (kappar + 1.) +
         0.5 * (sqr(wri[ivx]) + sqr(wri[ivy]) + sqr(wri[ivz]));

    rhobar = 0.5 * (wli[IDN] + wri[IDN]);
    cbar = sqrt(0.5 * (1. + (1. / kappar + 1. / kappal) / 2.) *
                (wli[IPR] + wri[IPR]) / rhobar);
    pbar = 0.5 * (wli[IPR] + wri[IPR]) +
           0.5 * (rhobar * cbar) * (wli[ivx] - wri[ivx]);
    ubar = 0.5 * (wli[ivx] + wri[ivx]) +
           0.5 / (rhobar * cbar) * (wli[IPR] - wri[IPR]);

    Real rdl = 1., rdr = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      rdl -= wli[n];
      rdr -= wri[n];
    }

    if (ubar > 0.) {
      flx(IDN, k, j, i) = ubar * wli[IDN] * rdl;
      for (int n = 1; n <= NVAPOR; ++n)
        flx(n, k, j, i) = ubar * wli[IDN] * wli[n];
      flx(ivx, k, j, i) = ubar * wli[IDN] * wli[ivx] + pbar;
      flx(ivy, k, j, i) = ubar * wli[IDN] * wli[ivy];
      flx(ivz, k, j, i) = ubar * wli[IDN] * wli[ivz];
      flx(IEN, k, j, i) = ubar * wli[IDN] * hl;
    } else {
      flx(IDN, k, j, i) = ubar * wri[IDN] * rdr;
      for (int n = 1; n <= NVAPOR; ++n)
        flx(n, k, j, i) = ubar * wri[IDN] * wri[n];
      flx(ivx, k, j, i) = ubar * wri[IDN] * wri[ivx] + pbar;
      flx(ivy, k, j, i) = ubar * wri[IDN] * wri[ivy];
      flx(ivz, k, j, i) = ubar * wri[IDN] * wri[ivz];
      flx(IEN, k, j, i) = ubar * wri[IDN] * hr;
    }

    // sedimentation flux
    for (int n = 0; n < NCLOUD; ++n) {
      Real vsed = ubar + pmicro->vsedf[dir](n, k, j, i);
      if (vsed > 0.) {
        pmicro->SetFluidMassFlux(dir, n, k, j, i, vsed * wli[IDN] * rdl);
      } else {
        pmicro->SetFluidMassFlux(dir, n, k, j, i, vsed * wri[IDN] * rdr);
      }
    }
  }
}
