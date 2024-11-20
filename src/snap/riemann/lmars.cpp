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

// exo3
#include <exo3/exo3.hpp>

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
  MeshBlock *pmb = pmy_block;

  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pthermo->GetGammad();
  Real wli[NHYDRO], wri[NHYDRO];

  for (int i = il; i <= iu; ++i) {
    // copy local variables
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }

    // correction for gamma
    Real kappal = 1. / (pthermo->GetGamma(wli) - 1.);
    Real kappar = 1. / (pthermo->GetGamma(wri) - 1.);

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

    // sedimentation flux
    // In Athena++ passive tracer module, the tracer mixing ratio is defined as
    // the ratio of the density of the tracer to the density of the "dry"
    // species. This ensures conservation of the tracer concentration even with
    // condensation However, the first component of the primitive variable is
    // the "total" density To get the denisty of the dry species, the "dry
    // mixing ratio" (rdl/rdr) is multiplied
    // FIXME: remove this?
    if (ivx == IVX) {
      if (i == iu) return;
      for (int n = 0; n < NCLOUD + NPRECIP; ++n) {
        Real rho = wri[IDN] * wri[1 + NVAPOR + n];
        auto vsed = pmicro->vsedf[0](n, k, j, i);
        flx(1 + NVAPOR + n, k, j, i) += vsed * rho;
        flx(ivx, k, j, i) += vsed * rho * wri[ivx];
        flx(ivy, k, j, i) += vsed * rho * wri[ivy];
        flx(ivz, k, j, i) += vsed * rho * wri[ivz];
        flx(IEN, k, j, i) += vsed * rho * hr;
      }
    }
  }
}
