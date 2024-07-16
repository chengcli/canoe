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

// exo3
#include <exo3/cubed_sphere.hpp>

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
  auto pcoord = pmy_block->pcoord;

  auto pthermo = Thermodynamics::GetInstance();
  auto pmicro = pmy_block->pimpl->pmicro;

  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pmy_block->peos->GetGamma();
  Real wli[NHYDRO], wri[NHYDRO];

  AthenaArray<Real> empty{};
#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of projection
  if (ivx == IVY) {
    pcoord->PrimToLocal2(k, j, il, iu, empty, wl, wr, empty);
  } else if (ivx == IVZ) {
    pcoord->PrimToLocal3(k, j, il, iu, empty, wl, wr, empty);
  }
#endif  // AFFINE or CUBED_SPHERE

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
    Real vy, vz, KE_l, KE_r;
#ifdef CUBED_SPHERE
    auto pexo3 = pmy_block->pimpl->pexo3;
    if (ivx == IVX) {
      pexo3->ContravariantVectorToCovariant(j, k, wli[IVY], wli[IVZ], &vy, &vz);
      KE_l = 0.5 * (SQR(wli[IVX]) + wli[IVY] * vy + wli[IVZ] * vz);

      pexo3->ContravariantVectorToCovariant(j, k, wri[IVY], wri[IVZ], &vy, &vz);
      KE_r = 0.5 * (SQR(wri[IVX]) + wri[IVY] * vy + wri[IVZ] * vz);
    } else {
      KE_l = 0.5 * (SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      KE_r = 0.5 * (SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
    }
#else   // NOT CUBED_SPHERE
    KE_l = 0.5 * (SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
    KE_r = 0.5 * (SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
#endif  // CUBED_SPHERE

    hl = wli[IPR] / wli[IDN] * (kappal + 1.) + KE_l;
    hr = wri[IPR] / wri[IDN] * (kappar + 1.) + KE_r;

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



#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of deprojection
  if (ivx == IVY) {
    pcoord->FluxToGlobal2(k, j, il, iu, empty, empty, flx, empty, empty);
  } else if (ivx == IVZ) {
    pcoord->FluxToGlobal3(k, j, il, iu, empty, empty, flx, empty, empty);
  }
#endif  // AFFINE or CUBED_SPHERE

#ifdef CUBED_SPHERE  // projection from contravariant fluxes to covariant
  Real x, y;

  if (ivx == IVX) {
    x = tan(pcoord->x2v(j));
    y = tan(pcoord->x3v(k));
  } else if (ivx == IVY) {
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
    Real ty = flx(IVY, k, j, i);
    Real tz = flx(IVZ, k, j, i);

    // Transform fluxes
    flx(IVY, k, j, i) = ty + tz * cth;
    flx(IVZ, k, j, i) = tz + ty * cth;
  }
#endif  // CUBED_SPHERE

    // sedimentation flux
    // In Athena++ passive tracer module, the tracer mixing ratio is defined as
    // the ratio of the density of the tracer to the density of the "dry"
    // species. This ensures conservation of the tracer concentration even with
    // condensation However, the first component of the primitive variable is
    // the "total" density To get the denisty of the dry species, the "dry
    // mixing ratio" (rdl/rdr) is multiplied
    for (int n = 0; n < NCLOUD; ++n) {
      Real vfld = ubar + pmicro->vsedf[dir](n, k, j, i);
      if (vfld > 0.) {
        pmicro->mass_flux[dir](n, k, j, i) = vfld * wli[IDN] * rdl;
      } else {
        pmicro->mass_flux[dir](n, k, j, i) = vfld * wri[IDN] * rdr;
      }
    }
  }
}
