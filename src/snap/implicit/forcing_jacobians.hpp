#ifndef SRC_SNAP_IMPLICIT_FORCING_JACOBIANS_HPP_
#define SRC_SNAP_IMPLICIT_FORCING_JACOBIANS_HPP_

// Athena++ header
#include <coordinates/coordinates.hpp>
#include <mesh/mesh.hpp>

// canoe headers
#include "implicit_solver.hpp"

inline double cot(double x) { return (1. / tan(x)); }

template <typename T>
void ImplicitSolver::JacobianGravityCoriolis(T &jac, Real const prim[], int k,
                                             int j, int i) {
  Real omega1 = 0, omega2 = 0, omega3 = 0, theta, phi;
  Real omegax = pmy_hydro->hsrc.GetOmegaX();
  Real omegay = pmy_hydro->hsrc.GetOmegaY();
  Real omegaz = pmy_hydro->hsrc.GetOmegaZ();
  Real grav = pmy_hydro->hsrc.GetG1();
  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;

  MeshBlock *pmb = pmy_hydro->pmy_block;

  jac.setZero();
  if (strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    omega1 = omegaz;
    omega2 = omegax;
    omega3 = omegay;
  } else if (strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    theta = pmb->pcoord->x2v(j);

    omega1 = cos(theta) * omegax + sin(theta) * omegay;
    omega2 = -sin(theta) * omegax + cos(theta) * omegay;
    omega3 = omegaz;
  } else if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    theta = pmb->pcoord->x2v(j);
    phi = pmb->pcoord->x3v(k);

    omega1 = sin(theta) * cos(phi) * omegax + sin(theta) * sin(phi) * omegay +
             cos(theta) * omegaz;
    omega2 = cos(theta) * cos(phi) * omegax + cos(theta) * sin(phi) * omegay -
             sin(theta) * omegaz;
    omega3 = -sin(phi) * omegax + cos(phi) * omegay;

    Real v1 = prim[IVX];
    Real v2 = prim[IVY];
    Real v3 = prim[IVZ];
    Real ss = v1 * v1 + v2 * v2 + v3 * v3;
    Real gm1 = pmb->peos->getGamma() - 1;
    Real r = pmb->pcoord->x1v(i);

    jac(ivx, idn) += (v1 * v1 + (gm1 - 1) * ss) / r;
    jac(ivx, ivx) += -2 * gm1 * v1 / r;
    jac(ivx, ivy) += 2 * (1 - gm1) * v2 / r;
    jac(ivx, ivz) += 2 * (1 - gm1) * v3 / r;
    jac(ivx, ien) += 2 * gm1 / r;

    jac(ivy, idn) += (v1 * v2 + (gm1 / 2. * ss - v3 * v3) * cot(theta)) / r;
    jac(ivy, ivx) += -(v2 + gm1 * v1 * cot(theta)) / r;
    jac(ivy, ivy) += -(v1 + gm1 * v2 * cot(theta)) / r;
    jac(ivy, ivz) += (2 - gm1) * v3 * cot(theta) / r;
    jac(ivy, ien) += gm1 * cot(theta) / r;

    jac(ivz, idn) += v3 * (v1 + v2 * cot(theta)) / r;
    jac(ivz, ivx) += -v3 / r;
    jac(ivz, ivy) += -v3 * cot(theta) / r;
    jac(ivz, ivz) += -(v1 + v2 * cot(theta)) / r;
  }

  jac(ivx, idn) += grav;
  jac(ien, ivx) += grav;
  jac(ivx, ivy) += 2. * omega3;
  jac(ivx, ivz) += -2 * omega2;
  jac(ivy, ivx) += -2 * omega3;
  jac(ivy, ivz) += 2 * omega1;
  jac(ivz, ivx) += 2 * omega2;
  jac(ivz, ivy) += -2 * omega1;

  if (mydir_ == X2DIR)
    jac = p2_ * jac * p3_;
  else if (mydir_ == X3DIR)
    jac = p3_ * jac * p2_;

  /*} else if (dir == X2DIR) {
    jac(ivz,idn) = grav;
    jac(ien,ivz) = grav;
    Real tmp = omega1;
    omega1 = omega2;
    omega2 = omega3;
    omega3 = tmp;
  } else { // X3DIR
    jac(ivy,idn) = grav;
    jac(ien,ivy) = grav;
    Real tmp = omega3;
    omega3 = omega2;
    omega2 = omega1;
    omega1 = tmp;
  }*/
}

#endif  // SRC_SNAP_IMPLICIT_FORCING_JACOBIANS_HPP_
