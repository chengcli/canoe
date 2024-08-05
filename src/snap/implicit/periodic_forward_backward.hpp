#ifndef SRC_SNAP_IMPLICIT_PERIODIC_FORWARD_BACKWARD_HPP_
#define SRC_SNAP_IMPLICIT_PERIODIC_FORWARD_BACKWARD_HPP_

// C/C++ headers
#include <vector>

// Eigen headers
#include <Eigen/Core>
#include <Eigen/Dense>

// Athena++ headers
#include "communication.hpp"

template <typename T1, typename T2>
void ImplicitSolver::PeriodicForwardSweep(std::vector<T1> &diag,
                                          std::vector<T1> &diagL,
                                          std::vector<T1> &diagU,
                                          std::vector<T2> &corr, Real dt, int k,
                                          int j, int il, int iu) {
  T1 sum_beta_gamma, phi;
  T2 rhs, sum_beta_zeta;
  std::vector<T1> beta(diag.size()), gamma(diag.size());
  std::vector<T1> ainv(diag.size());
  std::vector<T2> zeta(diag.size());

  // start periodic linear system solver, ref: El-Mikkawy (2005).
  // step 1: compute alpha, beta, gamma, delta from is to ie-1.
  //                                                 -1
  // alpha[1] = a[1], gamma[1] = U, beta[1] = L*alpha[1]

  // diag -> d
  // diagU -> a
  // diagL -> b
  // alpha -> c
  // gamma -> v
  // beta -> h
  // zeta -> k

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

  // if (last_block)
  //   SendBuffer(diagU[iu], k, j, tblock);
  // if (first_block)
  //   RecvBuffer(diagU[il-1], k, j, bblock);

  if (first_block) {
    // c[1] = d[1]
    ainv[il] = diag[il].inverse().eval();
    // v[1] = t (upper corner)
    gamma[il] = diagL[il];
    //             -1
    // h[1] = s*c[1]  (lower corner)
    beta[il] = diagU[il - 1] * ainv[il];
    // r[1] = k[1]
    zeta[il] = rhs;
    sum_beta_gamma.setZero();
    sum_beta_zeta.setZero();
  } else {
    RecvBuffer(ainv[il - 1], gamma[il - 1], beta[il - 1], zeta[il - 1],
               sum_beta_gamma, sum_beta_zeta, k, j, bblock);
    //                             -1
    // c[il] = d[il] - b[il]*c[il-1]*a[il-1]
    ainv[il] = diag[il] - diagL[il] * ainv[il - 1] * diagU[il - 1];
    ainv[il] = ainv[il].inverse().eval();
    //                      -1
    // v[il] = -b[il]*c[il-1]*v[il-1]
    gamma[il] = -diagL[il] * ainv[il - 1] * gamma[il - 1];
    //                              -1
    // h[il] = -h[il-1]*a[il-1]*c[il]
    beta[il] = -beta[il - 1] * diagU[il - 1] * ainv[il];
    //                           -1
    // r[il] = k[il] - b[il]*c[il]*r[il-1]
    zeta[il] = rhs - diagL[il] * ainv[il - 1] * zeta[il - 1];
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

    //                         -1
    // c[i] = d[i] - b[i]*c[i-1]*a[i-1]
    ainv[i] = diag[i] - diagL[i] * ainv[i - 1] * diagU[i - 1];
    ainv[i] = ainv[i].inverse().eval();
    //                   -1
    // v[i] = -b[i]*c[i-1]*v[i-1]
    gamma[i] = -diagL[i] * ainv[i - 1] * gamma[i - 1];
    //                          -1
    // h[i] = -h[i-1]*a[i-1]*c[i]
    beta[i] = -beta[i - 1] * diagU[i - 1] * ainv[i];
    //                       -1
    // r[i] = k[i] - b[i]*c[i]*r[i-1]
    zeta[i] = rhs - diagL[i] * ainv[i - 1] * zeta[i - 1];
  }

  if (last_block) {  // I'm the last block
    //                               -1
    // v[n-1] = a[n-1] - b[n-1]*c[n-2]*v[n-2];
    gamma[iu - 1] =
        diagU[iu - 1] - diagL[iu - 1] * ainv[iu - 2] * gamma[iu - 2];
    //                                      -1
    // h[n-1] = (b[n] - h[n-2]*a[n-2])*c[n-1]
    beta[iu - 1] = (diagL[iu] - beta[iu - 2] * diagU[iu - 2]) * ainv[iu - 1];
    for (int i = il; i < iu; ++i) {
      // h[i]*v[i]
      sum_beta_gamma += beta[i] * gamma[i];
      // h[i]*r[i]
      sum_beta_zeta += beta[i] * zeta[i];
    }
    //                 n-1
    // c[n] = d[n] - Sum {h[i]*v[i]}
    //                 i=1
    ainv[iu] = diag[iu] - sum_beta_gamma;
    ainv[iu] = ainv[iu].inverse().eval();
    //                 n-1
    // r[n] = k[n] - Sum {h[i]*r[i]}
    //                 i=1
    zeta[iu] = rhs - sum_beta_zeta;
  } else {
    for (int i = il; i <= iu; ++i) {
      // h[i]*v[i]
      sum_beta_gamma += beta[i] * gamma[i];
      // h[i]*r[i]
      sum_beta_zeta += beta[i] * zeta[i];
    }
    SendBuffer(ainv[iu], gamma[iu], beta[iu], zeta[iu], sum_beta_gamma,
               sum_beta_zeta, k, j, tblock);
  }

  SaveCoefficients(ainv, gamma, zeta, diagU, k, j, il, iu);
}

// delta -> x
// zeta -> r
// alpha -> c
// diagU -> a
// gamma -> v

template <typename T1, typename T2>
void ImplicitSolver::PeriodicBackwardSubstitution(std::vector<T1> &ainv,
                                                  std::vector<T1> &diagU,
                                                  std::vector<T2> &delta, int k,
                                                  int j, int il, int iu) {
  T2 delta_last;
  std::vector<T1> gamma(diagU.size());
  std::vector<T2> zeta(diagU.size());

  LoadCoefficients(ainv, gamma, zeta, diagU, k, j, il, iu);
  if (last_block) {  // I'm the last block
    // a[n-1] = 0
    diagU[iu - 1].setZero();
    //           -1
    // x[n] = c[n]*r[n]
    delta[iu] = ainv[iu] * zeta[iu];
    delta_last = delta[iu];
  } else {
    RecvBuffer(delta[iu + 1], delta_last, k, j, tblock);
    delta[iu] = ainv[iu] *
                (zeta[iu] - diagU[iu] * delta[iu + 1] - gamma[iu] * delta_last);
  }

  // backward substitution
  for (int i = iu - 1; i >= il; --i)
    //           -1
    // x[i] = c[i]*(r[i] - a[i]*x[i+1] - v[i]*x[n])
    delta[i] =
        ainv[i] * (zeta[i] - diagU[i] * delta[i + 1] - gamma[i] * delta_last);

  // update conserved variables
  for (int i = il; i <= iu; ++i) {
    if (T2::RowsAtCompileTime == 3) {  // partial matrix
      du_(IDN, k, j, i) = delta[i](0);
      du_(IVX + mydir_, k, j, i) = delta[i](1);
      du_(IEN, k, j, i) = delta[i](2);
    } else {  // full matrix
      du_(IDN, k, j, i) = delta[i](0);
      du_(IVX + mydir_, k, j, i) = delta[i](1);
      du_(IVX + (IVY - IVX + mydir_) % 3, k, j, i) = delta[i](2);
      du_(IVX + (IVZ - IVX + mydir_) % 3, k, j, i) = delta[i](3);
      du_(IEN, k, j, i) = delta[i](4);
    }
    for (int n = 1; n < IVX; ++n) du_(IDN, k, j, i) -= du_(n, k, j, i);
  }

  if (!first_block) SendBuffer(delta[il], delta_last, k, j, bblock);

#ifdef MPI_PARALLEL
  MPI_Status status;

  if (!last_block && (tblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_data6_[k][j], &status);

  if (!first_block && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_data2_[k][j], &status);
#endif
}

#endif  //  SRC_SNAP_IMPLICIT_PERIODIC_FORWARD_BACKWARD_HPP_
