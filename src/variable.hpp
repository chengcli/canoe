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
  enum { Size = NHYDRO + NSCALARS + NSTATIC };

  //! data pointers
  //! hydro data
  Real *w;

  //! cloud data
  Real *c;

  //! tracer data
  Real *s;

  //! chemistry data
  Real *q;

  //! static data
  Real *x;

  //! particle data
  Real *p;

  // constructor
  Variable() {
    w = data_.data();
    c = w + NHYDRO;
    s = c + NCLOUD;
    q = s + NTRACER;
    x = q + NCHEMISTRY;
  }

 private:
  // data holder
  std::array<Real, Size> data_;
};

#endif  // SRC_SNAP_VARIABLE_HPP_
