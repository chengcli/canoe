#ifndef SRC_TRANSPORT_STREAM_TRANSPORT_HPP_
#define SRC_TRANSPORT_STREAM_TRANSPORT_HPP_

// canoe
#include <configure.h>

// dealii headers
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

// forward declaration
template <typename>
class AthenaArray;

class StreamTransport {
 public:
  StreamTransport(int rows, int cols, bool fourth_order = true);
  ~StreamTransport();

  void setDiffusionMatrix(Real kdiff, Real dx);

  void setDiffusionMatrix(AthenaArray<Real> const& mkdiff);

  void setAdvectionMatrix(AthenaArray<Real> const& streamf, Real dx);

  void assembleSystem(Real dt, Real theta);

  // theta = 1: Backward Euler
  void evolve(AthenaArray<Real>* q, Real dt, Real theta = 1.);

  int64_t global(int i, int j) const { return i * cols_ + j; }

  int64_t globalh(int i, int j) const {
    return (i + GhostZoneSize) * colsh_ + (j + GhostZoneSize);
  }

 protected:
  // fourth order diffusion
  bool fourth_order_;
  // without ghost zones (halo)
  int rows_, cols_;
  int rank_;

  // with ghost zones (halo)
  int rowsh_, colsh_;
  int rankh_;

  // sparsity patterns
  dealii::SparsityPattern skh_, skk_, shk_;

  // diffusion matrix
  dealii::SparseMatrix<Real> diffusion_;

  // advection matrix
  dealii::SparseMatrix<Real> advection_;

  // stratch matrices, K-J, (K-J)*N
  dealii::SparseMatrix<Real> KmJ_, KmJmN_;

  // boundary conditions
  dealii::SparseMatrix<Real> bneumann_;

  // equation matrix and vectors
  dealii::Vector<Real> qvec_, dqvec_, rhs_;
  dealii::SparseMatrix<Real> mass_;

  // solver
  dealii::SolverControl control_;
  dealii::SolverGMRES<> solver_;
};

#endif  // SRC_TRANSPORT_STREAM_TRANSPORT_HPP_
