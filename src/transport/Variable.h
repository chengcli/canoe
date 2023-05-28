#ifndef VARIABLE_H_
#define VARIABLE_H_
#include <deal.II/lac/sparse_matrix.h>

#include <Eigen/Core>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Boundary.h"
#include "RectGrid.h"

template <typename _Scalar, int _Dim>
class Variable {};

template <typename _Scalar>
class Variable<_Scalar, 2>
    : public Eigen::Array<_Scalar, Eigen::Dynamic, Eigen::Dynamic> {
  template <typename STREAM>
  friend STREAM& operator<<(STREAM& os, const Variable& var) {
    os << std::setw(9) << var.m_name << "|";
    for (int i = 0; i < var.rows(); i++)
      os << std::setw(10) << i - var.m_grid.halo();
    os << std::endl;
    for (int i = -1; i < var.rows(); i++) os << "----------";
    os << std::endl;
    for (int j = 0; j < var.cols(); j++) {
      os << std::setw(9) << j - var.m_grid.halo();
      os << "|";
      for (int i = 0; i < var.rows(); i++) os << std::setw(10) << var(i, j);
      os << std::endl;
    }
    for (int i = -1; i < var.rows(); i++) os << "----------";
    os << std::endl;

    return os;
  }

 protected:
  typedef _Scalar Scalar;
  typedef Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> Base;
  enum { Dimension = 2 };

  const RectGrid<Scalar, Dimension>& m_grid;

  Eigen::Block<Base> m_value;

  std::string m_name, m_longname, m_units;

  std::array<BoundaryInfo<Scalar>, 4> m_boundary;

  dealii::SparsityPattern m_pattern;

  dealii::SparseMatrix<Scalar> m_neumann;

  dealii::Vector<Scalar> m_dirichlet;

 public:
  Variable() {}

  template <typename OtherDerived>
  Variable(const Eigen::DenseBase<OtherDerived>& other) : Base(other) {}

  template <typename OtherDerived>
  Variable& operator=(const Eigen::ArrayBase<OtherDerived>& other) {
    Base::operator=(other);
    return *this;
  }

  Variable(const Variable& other)
      : Base(other),
        m_grid(other.m_grid),
        m_value(*this, m_grid.halo(), m_grid.halo(), m_grid.rows(),
                m_grid.cols()),
        m_name(other.m_name),
        m_longname(other.m_longname),
        m_units(other.m_units),
        m_pattern(m_grid.rankh(), m_grid.rank(), 1),
        m_dirichlet(m_grid.rankh()) {
    for (int i = 0; i < m_grid.rows(); i++)
      for (int j = 0; j < m_grid.cols(); j++)
        m_pattern.add(m_grid.globalh(i, j), m_grid.global(i, j));

    for (int i = 0; i < m_grid.rankh(); i++) m_dirichlet(i) = 0.;
  }

  Variable& operator=(const Variable& other) {
    m_value = other.m_value;
    return *this;
  }

  Variable(const RectGrid<Scalar, Dimension>& grid, std::string name = "",
           std::string longname = "", std::string units = "")
      : Base(grid.rowsh(), grid.colsh()),
        m_grid(grid),
        m_value(*this, grid.halo(), grid.halo(), grid.rows(), grid.cols()),
        m_name(name),
        m_longname(longname),
        m_units(units),
        m_pattern(grid.rankh(), grid.rank(), 1),
        m_dirichlet(grid.rankh()) {
    for (int i = 0; i < m_grid.rows(); i++)
      for (int j = 0; j < m_grid.cols(); j++)
        m_pattern.add(m_grid.globalh(i, j), m_grid.global(i, j));

    for (int i = 0; i < m_grid.rankh(); i++) m_dirichlet(i) = 0.;
  }

  void set_boundary(int id, BoundaryType type, Scalar value) {
    if (id < 0 || id > 3) {
      std::cerr << "Boundary Id error" << std::endl;
      exit(1);
    }
    m_boundary[id].m_type = type;
    m_boundary[id].m_value = value;
  }

  void write_binary(const char* filename) const {
    std::ofstream out(filename,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    int nrows = m_grid.rows(), ncols = m_grid.cols(), ncolsh = m_grid.colsh(),
        nhalo = m_grid.halo();
    out.write((char*)(&nrows), sizeof(int));
    out.write((char*)(&ncols), sizeof(int));
    for (int i = nhalo; i < nhalo + nrows; i++)
      out.write((char*)(this->data() + nhalo + i * ncolsh),
                ncols * sizeof(Scalar));
    out.close();
  }

  void write_file(const char* filename) const {
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    for (int j = 0; j < m_grid.cols(); j++) {
      for (int i = 0; i < m_grid.rows(); i++)
        out << std::setw(14) << std::setprecision(6) << m_value(i, j);
      out << std::endl;
    }
  }

  void finish();

  const dealii::SparseMatrix<Scalar>& neumann() const { return m_neumann; }

  const dealii::Vector<Scalar>& dirichlet() const { return m_dirichlet; }

  inline Eigen::Block<Base>& val() { return m_value; }

  inline Scalar& val(int i, int j) { return m_value.coeffRef(i, j); }

  inline const Eigen::Block<Base>& val() const { return m_value; }

  inline const Scalar& val(int i, int j) const { return m_value.coeff(i, j); }

  inline std::string get_name() const { return m_name; }

  inline std::string get_longname() const { return m_longname; }

  inline std::string get_units() const { return m_units; }
};

template <typename _Scalar>
void Variable<_Scalar, 2>::finish() {
  switch (m_boundary[0].m_type) {
    case Dirichlet:
      for (int i = 0; i < m_grid.rows(); i++)
        m_dirichlet(m_grid.globalh(i, -1)) = m_boundary[0].m_value;
      break;
    case Neumann:
      break;
    case Periodic:
      break;
    default:
      std::cerr << "Boundary type not found" << std::endl;
      exit(1);
  }
  switch (m_boundary[1].m_type) {
    case Dirichlet:
      for (int i = 0; i < m_grid.rows(); i++)
        m_dirichlet(m_grid.globalh(i, m_grid.cols())) = m_boundary[1].m_value;
      break;
    case Neumann:
      break;
    case Periodic:
      break;
    default:
      std::cerr << "Boundary type not found" << std::endl;
      exit(1);
  }
  switch (m_boundary[2].m_type) {
    case Dirichlet:
      for (int j = 0; j < m_grid.cols(); j++)
        m_dirichlet(m_grid.globalh(-1, j)) = m_boundary[2].m_value;
      break;
    case Neumann:
      break;
    case Periodic:
      break;
    default:
      std::cerr << "Boundary type not found" << std::endl;
      exit(1);
  }
  switch (m_boundary[3].m_type) {
    case Dirichlet:
      for (int j = 0; j < m_grid.cols(); j++)
        m_dirichlet(m_grid.globalh(m_grid.rows(), j)) = m_boundary[3].m_value;
      break;
    case Neumann:
      break;
    case Periodic:
      break;
    default:
      std::cerr << "Boundary type not found" << std::endl;
      exit(1);
  }
  m_pattern.compress();

  m_neumann.reinit(m_pattern);
  for (int i = 0; i < m_grid.rows(); i++)
    for (int j = 0; j < m_grid.cols(); j++)
      m_neumann.set(m_grid.globalh(i, j), m_grid.global(i, j), 1.);
}

#endif
