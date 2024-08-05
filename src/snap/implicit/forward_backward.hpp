#ifndef SRC_SNAP_IMPLICIT_FORWARD_BACKWARD_HPP_
#define SRC_SNAP_IMPLICIT_FORWARD_BACKWARD_HPP_

// C/C++ headers
#include <vector>

// Eigen headers
#include <Eigen/Core>
#include <Eigen/Dense>

// Athena++ headers
#include "communication.hpp"

template <typename T1, typename T2>
void ImplicitSolver::ForwardSweep(std::vector<T1> &a, std::vector<T1> &b,
                                  std::vector<T1> &c, std::vector<T2> &delta,
                                  std::vector<T2> &corr, Real dt, int k, int j,
                                  int il, int iu) {
  T1 phi;
  T2 rhs;

  if (T2::RowsAtCompileTime == 3) {  // partial matrix
    rhs(0) = du_(IDN, k, j, il);
    for (int n = 1; n < IVX; ++n) rhs(0) += du_(n, k, j, il);
    rhs(0) /= dt;
    rhs(1) = du_(IVX + mydir_, k, j, il) / dt;
    rhs(2) = du_(IEN, k, j, il) / dt;
    rhs -= corr[il];
  } else {  // full matrix
    rhs(0) = du_(IDN, k, j, il);
    for (int n = 1; n < IVX; ++n) rhs(0) += du_(n, k, j, il);
    rhs(0) /= dt;
    rhs(1) = du_(IVX + mydir_, k, j, il) / dt;
    rhs(2) = du_(IVX + (IVY - IVX + mydir_) % 3, k, j, il) / dt;
    rhs(3) = du_(IVX + (IVZ - IVX + mydir_) % 3, k, j, il) / dt;
    rhs(4) = du_(IEN, k, j, il) / dt;
    // if (pmy_hydro->implicit_done != nullptr) {
    //   pmy_hydro->implicit_done->LoadForcingJacobian(phi,k,j,il,mydir_);
    //   rhs -= dt*phi*rhs;
    // }
  }

  if (!first_block) {
    RecvBuffer(a[il - 1], delta[il - 1], k, j, bblock);
    a[il] = (a[il] - b[il] * a[il - 1]).inverse().eval();
    delta[il] = a[il] * (rhs - b[il] * delta[il - 1]);
    a[il] *= c[il];
  } else {
    a[il] = a[il].inverse().eval();
    delta[il] = a[il] * rhs;
    a[il] *= c[il];
  }

  for (int i = il + 1; i <= iu; ++i) {
    if (T2::RowsAtCompileTime == 3) {  // partial matrix
      rhs(0) = du_(IDN, k, j, i);
      for (int n = 1; n < IVX; ++n) rhs(0) += du_(n, k, j, i);
      rhs(0) /= dt;
      rhs(1) = du_(IVX + mydir_, k, j, i) / dt;
      rhs(2) = du_(IEN, k, j, i) / dt;
      rhs -= corr[i];
    } else {
      rhs(0) = du_(IDN, k, j, i);
      for (int n = 1; n < IVX; ++n) rhs(0) += du_(n, k, j, i);
      rhs(0) /= dt;
      rhs(1) = du_(IVX + mydir_, k, j, i) / dt;
      rhs(2) = du_(IVX + (IVY - IVX + mydir_) % 3, k, j, i) / dt;
      rhs(3) = du_(IVX + (IVZ - IVX + mydir_) % 3, k, j, i) / dt;
      rhs(4) = du_(IEN, k, j, i) / dt;
      // if (pmy_hydro->implicit_done != nullptr) {
      //   pmy_hydro->implicit_done->LoadForcingJacobian(phi,k,j,i,mydir_);
      //   rhs -= dt*phi*rhs;
      // }
    }

    a[i] = (a[i] - b[i] * a[i - 1]).inverse().eval();
    delta[i] = a[i] * (rhs - b[i] * delta[i - 1]);
    a[i] *= c[i];
  }

  SaveCoefficients(a, delta, k, j, il, iu);

  if (!last_block) SendBuffer(a[iu], delta[iu], k, j, tblock);
}

template <typename T1, typename T2>
void ImplicitSolver::BackwardSubstitution(std::vector<T1> &a,
                                          std::vector<T2> &delta, int k, int j,
                                          int il, int iu) {
  LoadCoefficients(a, delta, k, j, il, iu);
  if (!last_block) {
    RecvBuffer(delta[iu + 1], k, j, tblock);
    delta[iu] -= a[iu] * delta[iu + 1];
  }

  // update solutions, i=iu
  for (int i = iu - 1; i >= il; --i) delta[i] -= a[i] * delta[i + 1];

  // 7. update conserved variables, i = iu
  auto &w = pmy_block_->phydro->w;
  for (int i = il; i <= iu; ++i) {
    Real dens = du_(IDN, k, j, i);
    for (int n = 1; n < IVX; ++n) dens += du_(n, k, j, i);
    dens = delta[i](0) - dens;

    if (T2::RowsAtCompileTime == 3) {  // partial matrix
      du_(IDN, k, j, i) = delta[i](0);
      for (int n = 1; n < IVX; ++n) {
        du_(n, k, j, i) += dens * w(n, k, j, i);
      }
      du_(IVX + mydir_, k, j, i) = delta[i](1);
      du_(IEN, k, j, i) = delta[i](2);
    } else {  // full matrix
      du_(IDN, k, j, i) = delta[i](0);
      for (int n = 1; n < IVX; ++n) {
        du_(n, k, j, i) += dens * w(n, k, j, i);
      }
      du_(IVX + mydir_, k, j, i) = delta[i](1);
      du_(IVX + (IVY - IVX + mydir_) % 3, k, j, i) = delta[i](2);
      du_(IVX + (IVZ - IVX + mydir_) % 3, k, j, i) = delta[i](3);
      du_(IEN, k, j, i) = delta[i](4);
    }
    for (int n = 1; n < IVX; ++n) du_(IDN, k, j, i) -= du_(n, k, j, i);
  }

  if (!first_block) SendBuffer(delta[il], k, j, bblock);

#ifdef MPI_PARALLEL
  MPI_Status status;

  if (!last_block && (tblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_data2_[k][j], &status);

  if (!first_block && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_data1_[k][j], &status);
#endif
}

#endif  // SRC_SNAP_IMPLICIT_FORWARD_BACKWARD_HPP_
