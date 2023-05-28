#ifndef RECTGRID_H_
#define RECTGRID_H_
#include <deal.II/lac/sparse_matrix.h>

#include <Eigen/Core>

template <typename _Scalar, int _Dim>
class RectGrid {};

template <typename _Scalar>
class RectGrid<_Scalar, 2> {
 protected:
  typedef _Scalar Scalar;
  enum { Dimension = 2 };

  int m_nrows, m_ncols, m_nhalo, m_nrowsh, m_ncolsh;

  long m_rank, m_rankh;

  Scalar m_dd;  // reference grid size

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> m_row_map, m_row_axis;

  Eigen::Matrix<Scalar, 1, Eigen::Dynamic> m_col_map, m_col_axis;

  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m_rowcol_map;

  dealii::SparsityPattern m_pattern;

 public:
  RectGrid(int nrows, int ncols, int nhalo)
      : m_nrows(nrows),
        m_ncols(ncols),
        m_nhalo(nhalo),
        m_nrowsh(nrows + 2 * nhalo),
        m_ncolsh(ncols + 2 * nhalo),
        m_rank(nrows * ncols),
        m_rankh(m_nrowsh * m_ncolsh),
        m_row_map(nrows),
        m_row_axis(nrows),
        m_col_map(ncols),
        m_col_axis(ncols),
        m_rowcol_map(nrows, ncols),
        m_pattern(m_rank, m_rankh, 9) {
    for (int i = 0; i < m_nrows; i++)
      for (int j = 0; j < m_ncols; j++) {
        long k = global(i, j);
        m_pattern.add(k, globalh(i, j));
        m_pattern.add(k, globalh(i - 1, j));
        m_pattern.add(k, globalh(i, j - 1));
        m_pattern.add(k, globalh(i + 1, j));
        m_pattern.add(k, globalh(i, j + 1));
        m_pattern.add(k, globalh(i - 1, j - 1));
        m_pattern.add(k, globalh(i - 1, j + 1));
        m_pattern.add(k, globalh(i + 1, j - 1));
        m_pattern.add(k, globalh(i + 1, j + 1));
      }
    m_pattern.compress();

    // reference RectGrid is defined at [0, 1] x [0, 1]
    m_dd = 1. / (nrows + 1);
  }

  inline const dealii::SparsityPattern& pattern() const { return m_pattern; }

  inline long global(int i, int j) const { return i * m_ncols + j; }

  inline long globalh(int i, int j) const {
    return (i + m_nhalo) * m_ncolsh + (j + m_nhalo);
  }

  inline int rows() const { return m_nrows; }

  inline int rowsh() const { return m_nrowsh; }

  inline int cols() const { return m_ncols; }

  inline int colsh() const { return m_ncolsh; }

  inline int halo() const { return m_nhalo; }

  inline long rank() const { return m_rank; }

  inline long rankh() const { return m_rankh; }

  inline Scalar dd() const { return m_dd; }
};

#endif
