#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include <deal.II/lac/sparse_matrix.h>

template <typename _Scalar, int _Dim, int _Order = 4>
class Diffusion {};

template <typename _Scalar>
class Diffusion<_Scalar, 2, 2> {
 protected:
  typedef _Scalar Scalar;
  enum { Dimension = 2 };

  const RectGrid<Scalar, Dimension>& m_grid;

  const Variable<Scalar, Dimension>& m_eddy;

  dealii::SparseMatrix<Scalar> m_jacobian;

 public:
  Diffusion(const RectGrid<Scalar, Dimension>& grid,
            const Variable<Scalar, Dimension>& eddy)
      : m_grid(grid), m_eddy(eddy), m_jacobian(grid.pattern()) {}

  void update() {
    for (int i = 0; i < m_grid.rows(); i++)
      for (int j = 0; j < m_grid.cols(); j++) {
        long k = m_grid.global(i, j);
        m_jacobian.set(k, m_grid.globalh(i, j), -4.);
        m_jacobian.set(k, m_grid.globalh(i - 1, j), 1.);
        m_jacobian.set(k, m_grid.globalh(i + 1, j), 1.);
        m_jacobian.set(k, m_grid.globalh(i, j - 1), 1.);
        m_jacobian.set(k, m_grid.globalh(i, j + 1), 1.);
      }
    m_jacobian *= m_eddy.val(0, 0) / (m_grid.dd() * m_grid.dd());
  }

  inline dealii::SparseMatrix<Scalar>& operator()() { return m_jacobian; }

  inline const dealii::SparseMatrix<Scalar>& operator()() const {
    return m_jacobian;
  }
};

template <typename _Scalar>
class Diffusion<_Scalar, 2, 4> {
 protected:
  typedef _Scalar Scalar;
  enum { Dimension = 2 };

  const RectGrid<Scalar, Dimension>& m_grid;

  const Variable<Scalar, Dimension>& m_eddy;

  dealii::SparseMatrix<Scalar> m_jacobian;

 public:
  Diffusion(const RectGrid<Scalar, Dimension>& grid,
            const Variable<Scalar, Dimension>& eddy)
      : m_grid(grid), m_eddy(eddy), m_jacobian(grid.pattern()) {}

  void update() {
    for (int i = 0; i < m_grid.rows(); i++)
      for (int j = 0; j < m_grid.cols(); j++) {
        long k = m_grid.global(i, j);
        m_jacobian.set(k, m_grid.globalh(i, j), -10. / 3.);
        m_jacobian.set(k, m_grid.globalh(i - 1, j), 2. / 3.);
        m_jacobian.set(k, m_grid.globalh(i, j - 1), 2. / 3.);
        m_jacobian.set(k, m_grid.globalh(i + 1, j), 2. / 3.);
        m_jacobian.set(k, m_grid.globalh(i, j + 1), 2. / 3.);
        m_jacobian.set(k, m_grid.globalh(i - 1, j - 1), 1. / 6.);
        m_jacobian.set(k, m_grid.globalh(i - 1, j + 1), 1. / 6.);
        m_jacobian.set(k, m_grid.globalh(i + 1, j - 1), 1. / 6.);
        m_jacobian.set(k, m_grid.globalh(i + 1, j + 1), 1. / 6.);
      }
    m_jacobian *= m_eddy.val(0, 0) / (m_grid.dd() * m_grid.dd());
  }

  inline dealii::SparseMatrix<Scalar>& operator()() { return m_jacobian; }

  inline const dealii::SparseMatrix<Scalar>& operator()() const {
    return m_jacobian;
  }
};

#endif
