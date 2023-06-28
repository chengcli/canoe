// C/C++
#include <iostream>
#include <vector>

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

// canoe
#include <impl.hpp>

// application
#include <application/application.hpp>

// athena
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include "../thermodynamics/thermodynamics.hpp"
#include "flux_decomposition.hpp"
#include "forward_backward.hpp"
#include "periodic_forward_backward.hpp"

void ImplicitSolver::PartialCorrection(AthenaArray<Real>& du,
                                       AthenaArray<Real> const& w, Real dt) {
  Application::Logger app("snap");
  app->Log("ImplicitSolver::PartialCorrection-X" + std::to_string(mydir_ + 1));

  int is, ie, js, je, ks, ke;
  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;

  MeshBlock const* pmb = pmy_block_;

  if (mydir_ == X1DIR) {
    ks = pmb->ks, js = pmb->js, is = pmb->is;
    ke = pmb->ke, je = pmb->je, ie = pmb->ie;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, k, j, i);
  } else if (mydir_ == X2DIR) {
    ks = pmb->is, js = pmb->ks, is = pmb->js;
    ke = pmb->ie, je = pmb->ke, ie = pmb->je;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, j, i, k);
  } else {  // X3DIR
    ks = pmb->js, js = pmb->is, is = pmb->ks;
    ke = pmb->je, je = pmb->ie, ie = pmb->ke;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, i, k, j);
  }

  // eigenvectors, eigenvalues, inverse matrix of eigenvectors.
  Eigen::Matrix<Real, 5, 5> Rmat, Lambda, Rimat;

  // reduced diffusion matrix |A_{i-1/2}|, |A_{i+1/2}|
  Eigen::Matrix<Real, 5, 5> Am, Ap, dfdq;
  Eigen::Matrix<Real, 3, 2> Am1, Ap1;
  Eigen::Matrix<Real, 3, 3> Am2, Ap2;
  Eigen::Matrix<Real, 3, 1> sm, sp;
  Eigen::Matrix<Real, 2, 1> dqm, dqp;

  Real prim[NHYDRO];  // Roe averaged primitive variables of cell i-1/2

  int nc = ie - is + 1 + 2 * NGHOST;
  std::vector<Eigen::Matrix<Real, 3, 3>> a(nc), b(nc), c(nc);
  std::vector<Eigen::Matrix<Real, 3, 1>> delta(nc);
  std::vector<Eigen::Matrix<Real, 3, 2>> dfdq1(nc);
  std::vector<Eigen::Matrix<Real, 3, 3>> dfdq2(nc);
  std::vector<Eigen::Matrix<Real, 3, 1>> corr(nc);

  Real* gamma_m1 = new Real[nc];

  // 0. forcing and volume matrix
  FindNeighbors();
  // SynchronizeConserved(du_, ks, ke, js, je, is, ie);
  // WaitToFinishSync(ks, ke, js, je, is, ie);

  Real gamma = pmy_block_->peos->GetGamma();
  Real grav = pmy_block_->phydro->hsrc.GetG1();
  Eigen::Matrix<Real, 3, 3> Phi, Dt, Bnd;
  if (mydir_ == X1DIR) {
    Phi << 0., 0., 0., grav, 0., 0., 0., grav, 0.;
  } else {
    Phi << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  }

  Dt << 1. / dt, 0., 0., 0., 1. / dt, 0., 0., 0., 1. / dt;

  Bnd << 1., 0., 0., 0., -1., 0., 0., 0., 1.;

  Real wl[NHYDRO], wr[NHYDRO];
  auto pthermo = pmy_block_->pimpl->pthermo;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      // calculate and save flux Jacobian matrix
      for (int i = is - 2; i <= ie + 1; ++i) {
        Real fsig = 1., feps = 1.;
        CopyPrimitives(wl, wr, w, k, j, i, mydir_);
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += wr[n] * (pthermo->GetCvRatio(n) - 1.);
          feps += wr[n] * (1. / pthermo->GetMassRatio(n) - 1.);
        }

        gamma_m1[i] = (gamma - 1.) * feps / fsig;
        // FluxJacobian(dfdq1[i], dfdq2[i], gamma_m1[i], w, k, j, i);
        FluxJacobian(dfdq, gamma_m1[i], wr, mydir_);
        dfdq1[i] << dfdq(idn, ivy), dfdq(idn, ivz), dfdq(ivx, ivy),
            dfdq(ivx, ivz), dfdq(ien, ivy), dfdq(ien, ivz);
        dfdq2[i] << dfdq(idn, idn), dfdq(idn, ivx), dfdq(idn, ien),
            dfdq(ivx, idn), dfdq(ivx, ivx), dfdq(ivx, ien), dfdq(ien, idn),
            dfdq(ien, ivx), dfdq(ien, ien);
      }

      // set up diffusion matrix and tridiagonal coefficients
      // left edge
      CopyPrimitives(wl, wr, w, k, j, is - 1, mydir_);
      Real gm1 = 0.5 * (gamma_m1[is - 2] + gamma_m1[is - 1]);
      RoeAverage(prim, gm1, wl, wr);
      Real cs = pmy_block_->peos->SoundSpeed(prim);
      Eigenvalue(Lambda, prim[IVX + mydir_], cs);
      Eigenvector(Rmat, Rimat, prim, cs, gm1, mydir_);
      Am = Rmat * Lambda * Rimat;
      Am1 << Am(idn, ivy), Am(idn, ivz), Am(ivx, ivy), Am(ivx, ivz),
          Am(ien, ivy), Am(ien, ivz);
      Am2 << Am(idn, idn), Am(idn, ivx), Am(idn, ien), Am(ivx, idn),
          Am(ivx, ivx), Am(ivx, ien), Am(ien, idn), Am(ien, ivx), Am(ien, ien);

      for (int i = is - 1; i <= ie; ++i) {
        CopyPrimitives(wl, wr, w, k, j, i + 1, mydir_);
        gm1 = 0.5 * (gamma_m1[i] + gamma_m1[i + 1]);
        RoeAverage(prim, gm1, wl, wr);
        Real cs = pmy_block_->peos->SoundSpeed(prim);
        Eigenvalue(Lambda, prim[IVX + mydir_], cs);
        Eigenvector(Rmat, Rimat, prim, cs, gm1, mydir_);
        Ap = Rmat * Lambda * Rimat;
        Ap1 << Ap(idn, ivy), Ap(idn, ivz), Ap(ivx, ivy), Ap(ivx, ivz),
            Ap(ien, ivy), Ap(ien, ivz);
        Ap2 << Ap(idn, idn), Ap(idn, ivx), Ap(idn, ien), Ap(ivx, idn),
            Ap(ivx, ivx), Ap(ivx, ien), Ap(ien, idn), Ap(ien, ivx),
            Ap(ien, ien);

        // set up diagonals a, b, c.
        Real aleft, aright, vol;
        Coordinates* pcoord = pmy_block_->pcoord;
        if (mydir_ == X1DIR) {
          aleft = pcoord->GetFace1Area(k, j, i);
          aright = pcoord->GetFace1Area(k, j, i + 1);
          vol = pcoord->GetCellVolume(k, j, i);
        } else if (mydir_ == X2DIR) {
          aleft = pcoord->GetFace2Area(j, i, k);
          aright = pcoord->GetFace2Area(j, i + 1, k);
          vol = pcoord->GetCellVolume(j, i, k);
        } else {  // X3DIR
          aleft = pcoord->GetFace3Area(i, k, j);
          aright = pcoord->GetFace3Area(i + 1, k, j);
          vol = pcoord->GetCellVolume(i, k, j);
        }
        a[i] = (Am2 * aleft + Ap2 * aright + (aright - aleft) * dfdq2[i]) /
                   (2. * vol) +
               Dt - Phi;
        b[i] = -(Am2 + dfdq2[i - 1]) * aleft / (2. * vol);
        c[i] = -(Ap2 - dfdq2[i + 1]) * aright / (2. * vol);

        // flux correction
        // dqm << du_(IVX+(IVY-IVX+mydir_)%3,k,j,i  ),
        // du_(IVX+(IVZ-IVX+mydir_)%3,k,j,i  ); dqp <<
        // du_(IVX+(IVY-IVX+mydir_)%3,k,j,i+1),
        // du_(IVX+(IVZ-IVX+mydir_)%3,k,j,i+1); sm = 0.5*((dfdq1[i-1] + Am1)*dqm
        // + (dfdq1[i] - Am1)*dqp); sp = 0.5*((dfdq1[i] + Ap1)*dqm + (dfdq1[i+1]
        // - Ap1)*dqp); corr[i] = (sp*aright - sm*aleft)/vol;
        corr[i].setZero();

        // Shift one cell: i -> i+1
        Am1 = Ap1;
        Am2 = Ap2;
      }

      // 5. fix boundary condition
      if (first_block && !periodic_boundary) a[is] += b[is] * Bnd;
      if (last_block && !periodic_boundary) a[ie] += c[ie] * Bnd;

      // 6. solve tridiagonal system using LU decomposition
      if (periodic_boundary)
        PeriodicForwardSweep(a, b, c, corr, dt, k, j, is, ie);
      else
        ForwardSweep(a, b, c, delta, corr, dt, k, j, is, ie);
    }

  if (periodic_boundary)
    PeriodicBackwardSubstitution(a, c, delta, ks, ke, js, je, is, ie);
  else
    BackwardSubstitution(a, delta, ks, ke, js, je, is, ie);

  if (mydir_ == X1DIR) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, k, j, i) = du_(IDN, k, j, i);
          du(IVX, k, j, i) = du_(IVX, k, j, i);
          du(IEN, k, j, i) = du_(IEN, k, j, i);
        }
  } else if (mydir_ == X2DIR) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, j, i, k) = du_(IDN, k, j, i);
          du(IVY, j, i, k) = du_(IVY, k, j, i);
          du(IEN, j, i, k) = du_(IEN, k, j, i);
        }
  } else {  // X3DIR
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, i, k, j) = du_(IDN, k, j, i);
          du(IVZ, i, k, j) = du_(IVZ, k, j, i);
          du(IEN, i, k, j) = du_(IEN, k, j, i);
        }
  }

  // pdebug->CheckConservation("du", du, is, ie, js, je, ks, ke);

  delete[] gamma_m1;
}
