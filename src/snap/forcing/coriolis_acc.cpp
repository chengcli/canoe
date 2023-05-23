//!  \brief source terms due to constant coriolis acceleration

// c++ header
#include <sstream>
#include <stdexcept>

// climath header
extern "C" {
#include <core.h>
}

// debugger
#include <debugger.hpp>

// Athena++ headers
#include <athena.hpp>
#include <coordinates/coordinates.hpp>
#include <hydro/hydro.hpp>
#include <hydro/srcterms/hydro_srcterms.hpp>
#include <mesh/mesh.hpp>

//! \brief add source terms for constant coriolis acceleration in
//  axial direction to conserved variables
void Forcing::Coriolis123(const Real dt, const AthenaArray<Real> *flx,
                          const AthenaArray<Real> &prim,
                          AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  pdebug->Enter("HydroSourceTerms::Coriolis123");

  if (omega1_ != 0.0 || omega2_ != 0.0 || omega3_ != 0.0) {
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real m1 = prim(IDN, k, j, i) * prim(IM1, k, j, i),
               m2 = prim(IDN, k, j, i) * prim(IM2, k, j, i),
               m3 = prim(IDN, k, j, i) * prim(IM3, k, j, i);
          cons(IM1, k, j, i) += 2. * dt * (omega3_ * m2 - omega2_ * m3);
          cons(IM2, k, j, i) += 2. * dt * (omega1_ * m3 - omega3_ * m1);
          if (pmb->ke > pmb->ks)  // 3d
            cons(IM3, k, j, i) += 2. * dt * (omega2_ * m1 - omega1_ * m2);
        }
      }
    }
  }

  pdebug->Leave();
  return;
}

//! \brief add source terms for constant coriolis acceleration in
//  cartesian coordinate to conserved variables
void Forcing::CoriolisXYZ(const Real dt, const AthenaArray<Real> *flx,
                          const AthenaArray<Real> &prim,
                          AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  pdebug->Enter("HydroSourceTerms::CoriolisXYZ");

  Real omega1, omega2, omega3, theta, phi;

  if (omegax_ != 0.0 || omegay_ != 0.0 || omegaz_ != 0.0) {
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          if (strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            omega1 = omegaz_;
            omega2 = omegax_;
            omega3 = omegay_;
          } else if (strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            theta = pmb->pcoord->x2v(j);

            omega1 = cos(theta) * omegax_ + sin(theta) * omegay_;
            omega2 = -sin(theta) * omegax_ + cos(theta) * omegay_;
            omega3 = omegaz_;
          } else if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            theta = pmb->pcoord->x2v(j);
            phi = pmb->pcoord->x3v(k);

            omega1 = sin(theta) * cos(phi) * omegax_ +
                     sin(theta) * sin(phi) * omegay_ + cos(theta) * omegaz_;
            omega2 = cos(theta) * cos(phi) * omegax_ +
                     cos(theta) * sin(phi) * omegay_ - sin(theta) * omegaz_;
            omega3 = -sin(phi) * omegax_ + cos(phi) * omegay_;
          } else {
            std::stringstream msg;
            msg << "### FATAL ERROR in HydroSourceTerms::CoriolisXYZ"
                << std::endl
                << "Unrecognized coordinate system" << std::endl;
            ATHENA_ERROR(msg);
          }

          Real m1 = prim(IDN, k, j, i) * prim(IM1, k, j, i),
               m2 = prim(IDN, k, j, i) * prim(IM2, k, j, i),
               m3 = prim(IDN, k, j, i) * prim(IM3, k, j, i);
          cons(IM1, k, j, i) += 2. * dt * (omega3 * m2 - omega2 * m3);
          cons(IM2, k, j, i) += 2. * dt * (omega1 * m3 - omega3 * m1);
          cons(IM3, k, j, i) += 2. * dt * (omega2 * m1 - omega1 * m2);
        }
    }
  }
  pdebug->Leave();

  return;
}
