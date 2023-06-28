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
  enum { Size = NHYDRO + NSCALARS + NCLOUDS + NCHEMISTRY + NSTATIC };

  //! data pointers
  //! hydro data
  Real *w;

  //! tracer data
  Real *s;

  //! cloud data
  Real *c;

  //! chemistry data
  Real *q;

  //! static data
  Real *x;

  //! particle data
  Real *p;

  // constructor
  Variable() {
    w = data_.data();
    s = w + NHYDRO;
    c = s + NSCALARS;
    q = c + NCLOUDS;
    x = q + NCHEMISTRY;
  }

 private:
  // data holder
  std::array<Real, Size> data_;
};

#endif  // SRC_SNAP_VARIABLE_HPP_
