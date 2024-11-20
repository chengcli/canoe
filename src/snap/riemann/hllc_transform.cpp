// C/C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// microphysics
#include <microphysics/microphysics.hpp>

// checks
#include <checks.hpp>

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//! \brief The HLLC Riemann solver for adiabatic hydrodynamics (use HLLE for
//! isothermal)

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
                          const int ivx, AthenaArray<Real> &wl,
                          AthenaArray<Real> &wr, AthenaArray<Real> &flx,
                          const AthenaArray<Real> &dxw) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pmicro = pmy_block->pimpl->pmicro;

  int ivy = IVX + ((ivx - IVX) + 1) % 3;
  int ivz = IVX + ((ivx - IVX) + 2) % 3;
  int dir = ivx - IVX;
  auto pcoord = pmy_block->pcoord;

  Real wli[(NHYDRO)], wri[(NHYDRO)];
  Real flxi[(NHYDRO)], fl[(NHYDRO)], fr[(NHYDRO)];
  Real gamma;
  if (GENERAL_EOS) {
    gamma = std::nan("");
  } else {
    gamma = pmy_block->peos->GetGamma();
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0 / gm1;

  AthenaArray<Real> empty{};
#if defined(AFFINE) || defined(CUBED_SPHERE)  // need of projection
  // if (ivx == IVX) {
  //   pcoord->PrimToLocal1(k, j, il, iu, empty, wl, wr, empty);
  if (ivx == IVY) {
    pcoord->PrimToLocal2(k, j, il, iu, empty, wl, wr, empty);
  } else if (ivx == IVZ) {
    pcoord->PrimToLocal3(k, j, il, iu, empty, wl, wr, empty);
  }
#endif  // AFFINE or CUBED_SPHERE

#pragma omp simd private(wli, wri, flxi, fl, fr)
#pragma distribute_point
  for (int i = il; i <= iu; ++i) {
    //--- Step 1.  Load L/R states into local variables
    wli[IDN] = wl(IDN, i);
    wli[IVX] = wl(ivx, i);
    wli[IVY] = wl(ivy, i);
    wli[IVZ] = wl(ivz, i);
    wli[IPR] = wl(IPR, i);

    wri[IDN] = wr(IDN, i);
    wri[IVX] = wr(ivx, i);
    wri[IVY] = wr(ivy, i);
    wri[IVZ] = wr(ivz, i);
    wri[IPR] = wr(IPR, i);

    for (int n = 1; n < IVX; ++n) {
      wli[n] = wl(n, i);
      wri[n] = wr(n, i);
    }

    // correction for gamma
    // left
    Real gammal = pthermo->GetGamma(wli);
    Real kappal = 1. / (gammal - 1.);

    // right
    Real gammar = pthermo->GetGamma(wri);
    Real kappar = 1. / (gammar - 1.);

    //--- Step 2.  Compute middle state estimates with PVRS (Toro 10.5.2)

    Real al, ar, el, er, vy, vz, KE_l, KE_r;
    Real cl = pmy_block->peos->SoundSpeed(wli);
    Real cr = pmy_block->peos->SoundSpeed(wri);

#ifdef CUBED_SPHERE
    auto pexo3 = pmy_block->pimpl->pexo3;
    if (ivx == IVX) {
      pexo3->ContravariantVectorToCovariant(j, k, wli[IVY], wli[IVZ], &vy, &vz);
      KE_l = 0.5 * wli[IDN] * (SQR(wli[IVX]) + wli[IVY] * vy + wli[IVZ] * vz);

      pexo3->ContravariantVectorToCovariant(j, k, wri[IVY], wri[IVZ], &vy, &vz);
      KE_r = 0.5 * wri[IDN] * (SQR(wri[IVX]) + wri[IVY] * vy + wri[IVZ] * vz);
    } else {
      KE_l = 0.5 * wli[IDN] * (SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
      KE_r = 0.5 * wri[IDN] * (SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
    }
#else   // NOT CUBED_SPHERE
    KE_l = 0.5 * wli[IDN] * (SQR(wli[IVX]) + SQR(wli[IVY]) + SQR(wli[IVZ]));
    KE_r = 0.5 * wri[IDN] * (SQR(wri[IVX]) + SQR(wri[IVY]) + SQR(wri[IVZ]));
#endif  // CUBED_SPHERE

    if (GENERAL_EOS) {
      el = pmy_block->peos->EgasFromRhoP(wli[IDN], wli[IPR]) + KE_l;
      er = pmy_block->peos->EgasFromRhoP(wri[IDN], wri[IPR]) + KE_r;
    } else {
      el = wli[IPR] * kappal + KE_l;
      er = wri[IPR] * kappar + KE_r;
    }
    Real rhoa = .5 * (wli[IDN] + wri[IDN]);  // average density
    Real ca = .5 * (cl + cr);                // average sound speed
    Real pmid = .5 * (wli[IPR] + wri[IPR] + (wli[IVX] - wri[IVX]) * rhoa * ca);
    Real umid =
        .5 * (wli[IVX] + wri[IVX] + (wli[IPR] - wri[IPR]) / (rhoa * ca));
    Real rhol = wli[IDN] + (wli[IVX] - umid) * rhoa / ca;  // mid-left density
    Real rhor = wri[IDN] + (umid - wri[IVX]) * rhoa / ca;  // mid-right density

    //--- Step 3.  Compute sound speed in L,R

    Real ql, qr;
    if (GENERAL_EOS) {
      Real gl = pmy_block->peos->AsqFromRhoP(rhol, pmid) * rhol / pmid;
      Real gr = pmy_block->peos->AsqFromRhoP(rhor, pmid) * rhor / pmid;
      ql = (pmid <= wli[IPR])
               ? 1.0
               : std::sqrt(1.0 + (gl + 1) / (2 * gl) * (pmid / wli[IPR] - 1.0));
      qr = (pmid <= wri[IPR])
               ? 1.0
               : std::sqrt(1.0 + (gr + 1) / (2 * gr) * (pmid / wri[IPR] - 1.0));
    } else {
      ql = (pmid <= wli[IPR]) ? 1.0
                              : std::sqrt(1.0 + (gammal + 1) / (2 * gammal) *
                                                    (pmid / wli[IPR] - 1.0));
      qr = (pmid <= wri[IPR]) ? 1.0
                              : std::sqrt(1.0 + (gammar + 1) / (2 * gammar) *
                                                    (pmid / wri[IPR] - 1.0));
    }

    //--- Step 4.  Compute the max/min wave speeds based on L/R

    al = wli[IVX] - cl * ql;
    ar = wri[IVX] + cr * qr;

    Real bp = ar > 0.0 ? ar : (TINY_NUMBER);
    Real bm = al < 0.0 ? al : -(TINY_NUMBER);

    //--- Step 5. Compute the contact wave speed and pressure

    Real vxl = wli[IVX] - al;
    Real vxr = wri[IVX] - ar;

    Real tl = wli[IPR] + vxl * wli[IDN] * wli[IVX];
    Real tr = wri[IPR] + vxr * wri[IDN] * wri[IVX];

    Real ml = wli[IDN] * vxl;
    Real mr = -(wri[IDN] * vxr);

    // Determine the contact wave speed...
    Real am = (tl - tr) / (ml + mr);
    // ...and the pressure at the contact surface
    Real cp = (ml * tr + mr * tl) / (ml + mr);
    cp = cp > 0.0 ? cp : 0.0;

    // No loop-carried dependencies anywhere in this loop
    //    #pragma distribute_point
    //--- Step 6. Compute L/R fluxes along the line bm, bp

    vxl = wli[IVX] - bm;
    vxr = wri[IVX] - bp;

    Real rdl = 1., rdr = 1.;
    for (int n = 1; n < IVX; ++n) {
      rdl -= wli[n];
      rdr -= wri[n];
    }

    fl[IDN] = wli[IDN] * vxl * rdl;
    fr[IDN] = wri[IDN] * vxr * rdr;

    for (int n = 1; n < IVX; ++n) {
      fl[n] = wli[IDN] * wli[n] * vxl;
      fr[n] = wri[IDN] * wri[n] * vxr;
    }

    fl[IVX] = wli[IDN] * wli[IVX] * vxl + wli[IPR];
    fr[IVX] = wri[IDN] * wri[IVX] * vxr + wri[IPR];

    fl[IVY] = wli[IDN] * wli[IVY] * vxl;
    fr[IVY] = wri[IDN] * wri[IVY] * vxr;

    fl[IVZ] = wli[IDN] * wli[IVZ] * vxl;
    fr[IVZ] = wri[IDN] * wri[IVZ] * vxr;

    fl[IEN] = el * vxl + wli[IPR] * wli[IVX];
    fr[IEN] = er * vxr + wri[IPR] * wri[IVX];

    //--- Step 8. Compute flux weights or scales

    Real sl, sr, sm;
    if (am >= 0.0) {
      sl = am / (am - bm);
      sr = 0.0;
      sm = -bm / (am - bm);
    } else {
      sl = 0.0;
      sr = -am / (bp - am);
      sm = bp / (bp - am);
    }

    //--- Step 9. Compute the HLLC flux at interface, including weighted
    // contribution
    // of the flux along the contact

    for (int n = 0; n < IVX; ++n) flxi[n] = sl * fl[n] + sr * fr[n];
    flxi[IVX] = sl * fl[IVX] + sr * fr[IVX] + sm * cp;
    flxi[IVY] = sl * fl[IVY] + sr * fr[IVY];
    flxi[IVZ] = sl * fl[IVZ] + sr * fr[IVZ];
    flxi[IEN] = sl * fl[IEN] + sr * fr[IEN] + sm * cp * am;

    for (int n = 0; n < IVX; ++n) flx(n, k, j, i) = flxi[n];
    flx(ivx, k, j, i) = flxi[IVX];
    flx(ivy, k, j, i) = flxi[IVY];
    flx(ivz, k, j, i) = flxi[IVZ];
    flx(IEN, k, j, i) = flxi[IEN];

    // sedimentation flux
    if (ivx == IVX) {
      if (i == iu) return;
      Real hr = (wri[IPR] * (kappar + 1.) + KE_r) / wri[IDN];
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

  check_hydro_riemann_solver_flux(flx, ivx, k, j, il, iu);
}
