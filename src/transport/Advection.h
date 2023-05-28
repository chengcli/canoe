#ifndef SRC_TRANSPORT_ADVECTION_H_
#define SRC_TRANSPORT_ADVECTION_H_
#include <deal.II/lac/sparse_matrix.h>

#include "Variable.h"

template <typename _Scalar, int _Dim, int _Order = 4>
class Advection {};

template <typename _Scalar>
class Advection<_Scalar, 2, 4> {
 protected:
  typedef _Scalar Scalar;
  enum { Dimension = 2 };

  const RectGrid<Scalar, Dimension>& m_grid;

  const Variable<Scalar, Dimension>& m_streamf;

  dealii::SparseMatrix<Scalar> m_jacobian;

 public:
  Advection(const RectGrid<Scalar, Dimension>& grid,
            const Variable<Scalar, Dimension>& streamf)
      : m_grid(grid), m_streamf(streamf), m_jacobian(grid.pattern()) {}

  void update() {
    for (int i = 0; i < m_grid.rows(); i++)
      for (int j = 0; j < m_grid.cols(); j++) {
        int64_t k = m_grid.global(i, j);
        m_jacobian.set(k, m_grid.globalh(i - 1, j - 1),
                       -m_streamf.val(i, j - 1) + m_streamf.val(i - 1, j));
        m_jacobian.set(k, m_grid.globalh(i - 1, j),
                       -m_streamf.val(i - 1, j - 1) - m_streamf.val(i, j - 1) +
                           m_streamf.val(i - 1, j + 1) +
                           m_streamf.val(i, j + 1));
        m_jacobian.set(k, m_grid.globalh(i - 1, j + 1),
                       +m_streamf.val(i, j + 1) - m_streamf.val(i - 1, j));
        m_jacobian.set(k, m_grid.globalh(i, j - 1),
                       -m_streamf.val(i + 1, j - 1) - m_streamf.val(i + 1, j) +
                           m_streamf.val(i - 1, j - 1) +
                           m_streamf.val(i - 1, j));
        m_jacobian.set(k, m_grid.globalh(i, j + 1),
                       +m_streamf.val(i + 1, j) + m_streamf.val(i + 1, j + 1) -
                           m_streamf.val(i - 1, j) -
                           m_streamf.val(i - 1, j + 1));
        m_jacobian.set(k, m_grid.globalh(i + 1, j - 1),
                       -m_streamf.val(i + 1, j) + m_streamf.val(i, j - 1));
        m_jacobian.set(k, m_grid.globalh(i + 1, j),
                       +m_streamf.val(i, j - 1) + m_streamf.val(i + 1, j - 1) -
                           m_streamf.val(i, j + 1) -
                           m_streamf.val(i + 1, j + 1));
        m_jacobian.set(k, m_grid.globalh(i + 1, j + 1),
                       +m_streamf.val(i + 1, j) - m_streamf.val(i, j + 1));
      }
    m_jacobian *= 1. / (12. * m_grid.dd() * m_grid.dd());
  }

  inline dealii::SparseMatrix<Scalar>& operator()() { return m_jacobian; }

  inline const dealii::SparseMatrix<Scalar>& operator()() const {
    return m_jacobian;
  }
};

#endif  // SRC_TRANSPORT_ADVECTION_H_
