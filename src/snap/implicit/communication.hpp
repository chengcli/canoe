#ifndef SRC_SNAP_IMPLICIT_COMMUNICATION_HPP_
#define SRC_SNAP_IMPLICIT_COMMUNICATION_HPP_

// C/C++
#include <cstring>
#include <functional>
#include <string>
#include <vector>

// athena
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <impl.hpp>

// snap
#include "implicit_solver.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

template <typename T>
void ImplicitSolver::SendBuffer(T const &a, int k, int j, NeighborBlock nb) {
  int s1 = a.size();
  // size_t phy = k << 10 | j << 2 | 0;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x1";

  memcpy(buffer_[k][j], a.data(), s1 * sizeof(Real));

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_block_->gid, phy);
    MPI_Isend(buffer_[k][j], s1, MPI_ATHENA_REAL, nb.snb.rank, tag,
              MPI_COMM_WORLD, &req_send_data1_[k][j]);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->pimpl->phevi->buffer_[k][j], buffer_[k][j],
                s1 * sizeof(Real));
  }
}

template <typename T1, typename T2>
void ImplicitSolver::SendBuffer(T1 const &a, T2 const &b, int k, int j,
                                NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size();
  // size_t phy = k << 10 | j << 2 | 1;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x2";

  memcpy(buffer_[k][j], a.data(), s1 * sizeof(Real));
  memcpy(buffer_[k][j] + s1, b.data(), s2 * sizeof(Real));

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_block_->gid, phy);
    MPI_Isend(buffer_[k][j], s1 + s2, MPI_ATHENA_REAL, nb.snb.rank, tag,
              MPI_COMM_WORLD, &req_send_data2_[k][j]);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->pimpl->phevi->buffer_[k][j], buffer_[k][j],
                (s1 + s2) * sizeof(Real));
  }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
void ImplicitSolver::SendBuffer(T1 const &a, T2 const &b, T3 const &c,
                                T4 const &d, T5 const &e, T6 const &f, int k,
                                int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size();
  // size_t phy = k << 10 | j << 2 | 2;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x6";

  Real *it = buffer_[k][j];
  memcpy(it, a.data(), s1 * sizeof(Real));
  it += s1;
  memcpy(it, b.data(), s2 * sizeof(Real));
  it += s2;
  memcpy(it, c.data(), s3 * sizeof(Real));
  it += s3;
  memcpy(it, d.data(), s4 * sizeof(Real));
  it += s4;
  memcpy(it, e.data(), s5 * sizeof(Real));
  it += s5;
  memcpy(it, f.data(), s6 * sizeof(Real));

  int st = s1 + s2 + s3 + s4 + s5 + s6;

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_block_->gid, phy);
    MPI_Isend(buffer_[k][j], st, MPI_ATHENA_REAL, nb.snb.rank, tag,
              MPI_COMM_WORLD, &req_send_data6_[k][j]);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->pimpl->phevi->buffer_[k][j], buffer_[k][j],
                st * sizeof(Real));
  }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7>
void ImplicitSolver::SendBuffer(T1 const &a, T2 const &b, T3 const &c,
                                T4 const &d, T5 const &e, T6 const &f,
                                T7 const &g, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size(), s7 = g.size();
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x7";

  Real *it = buffer_[k][j];
  memcpy(it, a.data(), s1 * sizeof(Real));
  it += s1;
  memcpy(it, b.data(), s2 * sizeof(Real));
  it += s2;
  memcpy(it, c.data(), s3 * sizeof(Real));
  it += s3;
  memcpy(it, d.data(), s4 * sizeof(Real));
  it += s4;
  memcpy(it, e.data(), s5 * sizeof(Real));
  it += s5;
  memcpy(it, f.data(), s6 * sizeof(Real));
  it += s6;
  memcpy(it, g.data(), s7 * sizeof(Real));

  int st = s1 + s2 + s3 + s4 + s5 + s6 + s7;

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_block_->gid, phy);
    MPI_Isend(buffer_[k][j], st, MPI_ATHENA_REAL, nb.snb.rank, tag,
              MPI_COMM_WORLD, &req_send_data7_[k][j]);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_block_->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->pimpl->phevi->buffer_[k][j], buffer_[k][j],
                st * sizeof(Real));
  }
}

template <typename T>
void ImplicitSolver::RecvBuffer(T &a, int k, int j, NeighborBlock nb) {
  int s1 = a.size();
  // size_t phy = k << 10 | j << 2 | 0;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x1";
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_block_->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_[k][j], s1, MPI_ATHENA_REAL, nb.snb.rank, tag,
             MPI_COMM_WORLD, &status);
#endif
  }  // local boundary

  memcpy(a.data(), buffer_[k][j], s1 * sizeof(Real));
}

template <typename T1, typename T2>
void ImplicitSolver::RecvBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size();
  // size_t phy = k << 10 | j << 2 | 1;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x2";
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_block_->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_[k][j], s1 + s2, MPI_ATHENA_REAL, nb.snb.rank, tag,
             MPI_COMM_WORLD, &status);
#endif
  }  // local boundary

  memcpy(a.data(), buffer_[k][j], s1 * sizeof(Real));
  memcpy(b.data(), buffer_[k][j] + s1, s2 * sizeof(Real));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
void ImplicitSolver::RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, T5 &e, T6 &f, int k,
                                int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size();
  // size_t phy = k << 10 | j << 2 | 2;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x6";
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  int st = s1 + s2 + s3 + s4 + s5 + s6;

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_block_->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_[k][j], st, MPI_ATHENA_REAL, nb.snb.rank, tag,
             MPI_COMM_WORLD, &status);
#endif
  }  // local boundary

  Real *it = buffer_[k][j];
  memcpy(a.data(), it, s1 * sizeof(Real));
  it += s1;
  memcpy(b.data(), it, s2 * sizeof(Real));
  it += s2;
  memcpy(c.data(), it, s3 * sizeof(Real));
  it += s3;
  memcpy(d.data(), it, s4 * sizeof(Real));
  it += s4;
  memcpy(e.data(), it, s5 * sizeof(Real));
  it += s5;
  memcpy(f.data(), it, s6 * sizeof(Real));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7>
void ImplicitSolver::RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, T5 &e, T6 &f, T7 &g,
                                int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size(), s7 = g.size();
  // size_t phy = k << 10 | j << 2 | 3;
  std::string phy = std::to_string(k) + "x" + std::to_string(j) + "x7";
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  int st = s1 + s2 + s3 + s4 + s5 + s6 + s7;

  if (nb.snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_block_->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_[k][j], st, MPI_ATHENA_REAL, nb.snb.rank, tag,
             MPI_COMM_WORLD, &status);
#endif
  }  // local boundary

  Real *it = buffer_[k][j];
  memcpy(a.data(), it, s1 * sizeof(Real));
  it += s1;
  memcpy(b.data(), it, s2 * sizeof(Real));
  it += s2;
  memcpy(c.data(), it, s3 * sizeof(Real));
  it += s3;
  memcpy(d.data(), it, s4 * sizeof(Real));
  it += s4;
  memcpy(e.data(), it, s5 * sizeof(Real));
  it += s5;
  memcpy(f.data(), it, s6 * sizeof(Real));
  it += s6;
  memcpy(g.data(), it, s7 * sizeof(Real));
}

template <typename T1, typename T2>
void ImplicitSolver::SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                                      int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size();
    memcpy(&coefficients_(k, j, i, 0), a[i].data(), s1 * sizeof(Real));
    memcpy(&coefficients_(k, j, i, s1), b[i].data(), s2 * sizeof(Real));
  }
}

template <typename T1, typename T2, typename T3, typename T4>
void ImplicitSolver::SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                                      std::vector<T3> &c, std::vector<T4> &d,
                                      int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size(), s3 = c[i].size(), s4 = d[i].size();
    memcpy(&coefficients_(k, j, i, 0), a[i].data(), s1 * sizeof(Real));
    memcpy(&coefficients_(k, j, i, s1), b[i].data(), s2 * sizeof(Real));
    memcpy(&coefficients_(k, j, i, s1 + s2), c[i].data(), s3 * sizeof(Real));
    memcpy(&coefficients_(k, j, i, s1 + s2 + s3), d[i].data(),
           s4 * sizeof(Real));
  }
}

template <typename T1, typename T2>
void ImplicitSolver::LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                                      int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size();
    memcpy(a[i].data(), &coefficients_(k, j, i, 0), s1 * sizeof(Real));
    memcpy(b[i].data(), &coefficients_(k, j, i, s1), s2 * sizeof(Real));
  }
}

template <typename T1, typename T2, typename T3, typename T4>
void ImplicitSolver::LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
                                      std::vector<T3> &c, std::vector<T4> &d,
                                      int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size(), s3 = c[i].size(), s4 = d[i].size();
    memcpy(a[i].data(), &coefficients_(k, j, i, 0), s1 * sizeof(Real));
    memcpy(b[i].data(), &coefficients_(k, j, i, s1), s2 * sizeof(Real));
    memcpy(c[i].data(), &coefficients_(k, j, i, s1 + s2), s3 * sizeof(Real));
    memcpy(d[i].data(), &coefficients_(k, j, i, s1 + s2 + s3),
           s4 * sizeof(Real));
  }
}

/*template<typename T>
void inline ImplicitSolver::SaveForcingJacobian(T &phi, int k, int j ,int i) {
  if (mydir == X1DIR) {
  } else if (mydir == X2DIR) {
    memcpy(jacobian_[j][i][k], phi.data(), phi.size()*sizeof(Real));
  } else {
    memcpy(jacobian_[i][k][j], phi.data(), phi.size()*sizeof(Real));
  }
}

template<typename T>
void inline ImplicitSolver::LoadForcingJacobian(T &phi, int k, int j ,int i,
  CoordinateDirection dir) {
  Eigen::Matrix<Real,5,5> tmp;

  if (dir == X1DIR) {
    memcpy(tmp.data(), jacobian_[k][j][i], tmp.size()*sizeof(Real));
  } else if (dir == X2DIR) {
    memcpy(tmp.data(), jacobian_[j][i][k], tmp.size()*sizeof(Real));
    tmp = p2_*tmp*p3_;
  } else { // X3DIR
    memcpy(tmp.data(), jacobian_[i][k][j], tmp.size()*sizeof(Real));
    tmp = p3_*tmp*p2_;
  }

  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      phi(i,j) = tmp(i,j);
}*/

#endif  // SRC_SNAP_IMPLICIT_COMMUNICATION_HPP_
