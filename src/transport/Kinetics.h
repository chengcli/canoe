#ifndef SRC_TRANSPORT_KINETICS_H_
#define SRC_TRANSPORT_KINETICS_H_

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "Advection.h"
#include "Diffusion.h"

class Kinetics {
  typedef double Scalar;
  enum { Dimension = 2 };

 protected:
  RectGrid<Scalar, Dimension> m_grid;

  Variable<Scalar, 2> m_psi, m_eddy, m_q;

  Advection<Scalar, 2, 4> m_advection;

  Diffusion<Scalar, 2, 4> m_diffusion;

  dealii::SparsityPattern m_pattern;

  dealii::SparseMatrix<Scalar> m_mass, m_force, m_adj;

  dealii::Vector<Scalar> m_src, m_rhs, m_bnd;

 public:
  Kinetics(int, int, int = 1);

  void initialize();

  void assemble(double, double);

  void checkout();

  void run(int);

 private:
  dealii::SparseMatrix<Scalar> m_buffer;
};

#endif  // SRC_TRANSPORT_KINETICS_H_
