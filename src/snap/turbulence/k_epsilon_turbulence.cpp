//! \file k_epsilon_turbulence.cpp
//  \brief implement K-Epsilon turbulence

// C/C++ headers
#include <algorithm>  // min,max

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include "turbulence_model.hpp"

KEpsilonTurbulence::KEpsilonTurbulence(MeshBlock *pmb, ParameterInput *pin)
    : TurbulenceModel(pmb, pin) {
  cmu_ = pin->GetOrAddReal("turbulence", "kepsilon.cmu", 0.09);
  c1_ = pin->GetOrAddReal("turbulence", "kepsilon.c1", 1.44);
  c2_ = pin->GetOrAddReal("turbulence", "kepsilon.c2", 1.92);
  sigk_ = pin->GetOrAddReal("turbulence", "kepsilon.sigk", 1.0);
  sige_ = pin->GetOrAddReal("turbulence", "kepsilon.sige", 1.3);

  // velocity scale, v ~ 1 m/s, tke ~ v^2 ~ 1 m^2/s^2
  Real tke0 = pin->GetOrAddReal("turbulence", "kepsilon.tke0", 1.);

  int il = pmb->is - NGHOST;
  int jl = pmb->js;
  int kl = pmb->ks;
  int iu = pmb->ie + NGHOST;
  int ju = pmb->je;
  int ku = pmb->ke;
  if (pmb->block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  if (pmb->block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        w(1, k, j, i) = tke0;
        // eps ~ tke/T ~ v^2/(dx/v) ~ v^3/dx
        Real volume = pmb->pcoord->GetCellVolume(k, j, i);
        Real eps0 = pow(tke0, 1.5) / pow(volume, 1. / 3.);
        w(0, k, j, i) = eps0;
      }
}

void KEpsilonTurbulence::Initialize() {
  auto pmb = pmy_block;
  auto phydro = pmb->phydro;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real rhod = phydro->w(IDN, k, j, i);
        u(0, k, j, i) = w(0, k, j, i) * rhod;
        u(1, k, j, i) = w(1, k, j, i) * rhod;
        mut(k, j, i) =
            cmu_ * rhod * w(1, k, j, i) * w(1, k, j, i) / w(0, k, j, i);
      }
}

inline Real Laplace_(AthenaArray<Real> const &mut, AthenaArray<Real> const &v,
                     int k, int j, int i, Coordinates *pcoord) {
  Real result = 0.;
  Mesh *pm = pcoord->pmy_block->pmy_mesh;

  // first dimension
  Real mut_p = (mut(k, j, i + 1) + mut(k, j, i)) / 2.;
  Real mut_m = (mut(k, j, i) + mut(k, j, i - 1)) / 2.;
  Real gradv_p =
      (v(k, j, i + 1) - v(k, j, i)) / (pcoord->x1v(i + 1) - pcoord->x1v(i));
  Real gradv_m =
      (v(k, j, i) - v(k, j, i - 1)) / (pcoord->x1v(i) - pcoord->x1v(i - 1));
  result += (mut_p * gradv_p - mut_m * gradv_m) / pcoord->dx1f(i);

  // second dimension
  if (pm->f2) {
    mut_p = (mut(k, j + 1, i) + mut(k, j, i)) / 2.;
    mut_m = (mut(k, j, i) + mut(k, j - 1, i)) / 2.;
    gradv_p =
        (v(k, j + 1, i) - v(k, j, i)) / (pcoord->x2v(j + 1) - pcoord->x2v(j));
    gradv_m =
        (v(k, j, i) - v(k, j - 1, i)) / (pcoord->x2v(j) - pcoord->x2v(j - 1));
    result += (mut_p * gradv_p - mut_m * gradv_m) / pcoord->dx2f(j);
  }

  // third dimension
  if (pm->f3) {
    mut_p = (mut(k + 1, j, i) + mut(k, j, i)) / 2.;
    mut_m = (mut(k, j, i) + mut(k - 1, j, i)) / 2.;
    gradv_p =
        (v(k + 1, j, i) - v(k, j, i)) / (pcoord->x3v(k + 1) - pcoord->x3v(k));
    gradv_m =
        (v(k, j, i) - v(k - 1, j, i)) / (pcoord->x3v(k) - pcoord->x3v(k - 1));
    result += (mut_p * gradv_p - mut_m * gradv_m) / pcoord->dx3f(k);
  }

  return result;
}

inline Real ShearProduction_(AthenaArray<Real> const &w, int k, int j, int i,
                             Coordinates *pcoord) {
  Real result = 0.;
  Mesh *pm = pcoord->pmy_block->pmy_mesh;
  Real gradv, gradv1, gradv2, gradv3, gradv4;

  // first dimension
  gradv = (w(IM1, k, j, i + 1) - w(IM1, k, j, i - 1)) /
          (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));
  result += 2 * gradv * gradv;

  // second dimension
  if (pm->f2) {
    gradv = (w(IM2, k, j + 1, i) - w(IM2, k, j - 1, i)) /
            (pcoord->x2v(j + 1) - pcoord->x2v(j - 1));
    gradv1 = (w(IM1, k, j + 1, i) - w(IM1, k, j - 1, i)) /
             (pcoord->x2v(j + 1) - pcoord->x2v(j - 1));
    gradv2 = (w(IM2, k, j, i + 1) - w(IM2, k, j, i - 1)) /
             (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));
    result += 2 * gradv * gradv + (gradv1 + gradv2) * (gradv1 + gradv2);
  }

  // thrid dimension
  if (pm->f3) {
    gradv = (w(IM3, k + 1, j, i) - w(IM3, k - 1, j, i)) /
            (pcoord->x3v(k + 1) - pcoord->x3v(k - 1));
    gradv1 = (w(IM1, k + 1, j, i) - w(IM1, k - 1, j, i)) /
             (pcoord->x3v(k + 1) - pcoord->x3v(k - 1));
    gradv2 = (w(IM3, k, j, i + 1) - w(IM3, k, j, i - 1)) /
             (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));
    gradv3 = (w(IM2, k + 1, j, i) - w(IM2, k - 1, j, i)) /
             (pcoord->x3v(k + 1) - pcoord->x3v(k - 1));
    gradv4 = (w(IM3, k, j + 1, i) - w(IM3, k, j - 1, i)) /
             (pcoord->x2v(j + 1) - pcoord->x2v(j - 1));
    result += 2 * gradv * gradv + (gradv1 + gradv2) * (gradv1 + gradv2) +
              (gradv3 + gradv4) * (gradv3 + gradv4);
  }

  return result;
}

void KEpsilonTurbulence::DriveTurbulence(Real dt) {
  // std::cout << "driving turbulence" << std::endl;
  auto pmb = pmy_block;
  auto phydro = pmb->phydro;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  AthenaArray<Real> eps, tke;
  eps.InitWithShallowSlice(w, 4, 0, 1);
  tke.InitWithShallowSlice(w, 4, 1, 1);

  Real s1, s2, s3;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real rhod = phydro->w(IDN, k, j, i);
        // shear production
        Real shear = ShearProduction_(w, k, j, i, pmb->pcoord);
        // std::cout << "shear = " << shear << std::endl;

        // turbulent dissipation, de/dt, eq2.2-1
        s1 = c1_ * mut(k, j, i) * eps(k, j, i) / tke(k, j, i) * shear;
        s2 = Laplace_(mut, eps, k, j, i, pmb->pcoord) / sige_;
        s3 = -c2_ * eps(k, j, i) * eps(k, j, i) / tke(k, j, i) * rhod;
        u(0, k, j, i) += (s1 + s2) * dt;
        if (u(0, k, j, i) + s3 * dt > 0)
          u(0, k, j, i) += s3 * dt;
        else
          u(0, k, j, i) /= 2.;
        // std::cout << "u = " << u(0,k,j,i) << std::endl;

        // turbulent kinetic energy, dk/dt, eq 2.2-2
        s1 = mut(k, j, i) * shear;
        s2 = Laplace_(mut, tke, k, j, i, pmb->pcoord) / sigk_;
        s3 = -eps(k, j, i) * rhod;
        u(1, k, j, i) += (s1 + s2) * dt;
        if (u(1, k, j, i) + s3 * dt > 0.)
          u(1, k, j, i) += s3 * dt;
        else
          u(1, k, j, i) /= 2.;
      }
}

void KEpsilonTurbulence::SetDiffusivity(AthenaArray<Real> &nu,
                                        AthenaArray<Real> &kappa,
                                        const AthenaArray<Real> &prim,
                                        const AthenaArray<Real> &bcc, int il,
                                        int iu, int jl, int ju, int kl,
                                        int ku) {
  AthenaArray<Real> eps, tke;
  eps.InitWithShallowSlice(w, 4, 0, 1);
  tke.InitWithShallowSlice(w, 4, 1, 1);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
#pragma omp simd
      for (int i = il; i <= iu; ++i) {
        // dynamic turbulent diffusivity
        mut(k, j, i) = cmu_ * prim(IDN, k, j, i) * tke(k, j, i) * tke(k, j, i) /
                       eps(k, j, i);

        // kinematic turbulent viscosity, mu_t = c_mu*k^2/epsilon, eq2.2-3
        nu(HydroDiffusion::DiffProcess::iso, k, j, i) =
            mut(k, j, i) / prim(IDN, k, j, i);
        kappa(HydroDiffusion::DiffProcess::iso, k, j, i) =
            nu(HydroDiffusion::DiffProcess::iso, k, j, i);
      }
    }
  }
}
