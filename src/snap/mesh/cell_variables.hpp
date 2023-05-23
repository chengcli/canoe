#ifndef CELL_VARIABLES_HPP
#define CELL_VARIABLES_HPP

// C/C++ header
#include <array>

// canoe header
#include <athena/athena.hpp>
#include <configure.hpp>

//! \class CellVariables
//  \brief a collection of all physical data in a computational cell
class CellVariables {
 public:
  enum { Size = NHYDRO + NSCALARS + NCLOUDS };

  //! data pointers
  //! hydro data
  Real *w;

  //! tracer data
  Real *s;

  //! cloud data
  Real *q;

  //! chemistry data
  Real *c;

  //! particle data
  Real *p;

  // constructor
  CellVariables() {
    w = data_.data();
    s = w + NHYDRO;
    q = s + NSCALARS;
    c = q + NCLOUDS;
  }

 private:
  // data holder
  std::array<Real, Size> data_;
};

#endif
