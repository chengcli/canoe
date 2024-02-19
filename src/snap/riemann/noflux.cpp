//! \file  noflux.cpp
//  \brief disable hydrodynamic flux

// Athena++
#include <athena/hydro/hydro.hpp>
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>

// canoe
#include <configure.hpp>
#include <exo3/cubed_sphere.hpp>
#include <impl.hpp>

// microphysics
#include <microphysics/microphysics.hpp>

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr,
                          AthenaArray<Real> &flx, const AthenaArray<Real> &dxw)
{
  int dir = ivx - IVX;
  Real pbar, rhobar, cbar;
  Real wli[NHYDRO], wri[NHYDRO];
  Real gamma = pmy_block->peos->GetGamma();

  AthenaArray<Real> empty{};
  auto pcoord = pmy_block->pcoord;

#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of projection
  // if (ivx == IVX) {
  //   pcoord->PrimToLocal1(k, j, il, iu, empty, wl, wr, empty);
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

    rhobar = 0.5 * (wli[IDN] + wri[IDN]);
    cbar = sqrt(0.5 * gamma *
                (wli[IPR] + wri[IPR]) / rhobar);
    pbar = 0.5 * (wli[IPR] + wri[IPR]) + 
           0.5 * (rhobar * cbar) * (wli[ivx] - wri[ivx]);

    for (int n = 0; n < NHYDRO; ++n)
      flx(n, k, j, i) = 0.;
    flx(ivx, k, j, i) = pbar;

    Real rdl = 1., rdr = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      rdl -= wli[n];
      rdr -= wri[n];
    }

    // sedimentation flux
    auto pmicro = pmy_block->pimpl->pmicro;
    for (int n = 0; n < NCLOUD; ++n) {
      Real vfld = pmicro->vsedf[dir](n, k, j, i);

      if (vfld > 0.) {
        pmicro->mass_flux[dir](n, k, j, i) = vfld * wli[IDN] * rdl;
      } else {
        pmicro->mass_flux[dir](n, k, j, i) = vfld * wri[IDN] * rdr;
      }
    }
  }

#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of deprojection
  // if (ivx == IVX) {
  //   pcoord->FluxToGlobal1(k, j, il, iu, empty, empty, flx, empty, empty);
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
}
