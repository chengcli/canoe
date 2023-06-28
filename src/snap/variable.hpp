#ifndef SRC_SNAP_VARIABLE_HPP_
#define SRC_SNAP_VARIABLE_HPP_

// C/C++
#include <array>

// canoe
#include <configure.hpp>

// athena
#include <athena/athena.hpp>

//! \class Variable
//  \brief a collection of all physical data in a computational cell
class Variable {
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
  Variable() {
    w = data_.data();
    s = w + NHYDRO;
    q = s + NSCALARS;
    c = q + NCLOUDS;
  }

 private:
  // data holder
  std::array<Real, Size> data_;
};

#endif  // SRC_SNAP_VARIABLE_HPP_
