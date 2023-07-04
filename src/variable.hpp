#ifndef SRC_VARIABLE_HPP_
#define SRC_VARIABLE_HPP_

// C/C++
#include <array>

// canoe
#include <configure.hpp>

// athena
#include <athena/athena.hpp>

//! \class Variable
//  \brief a collection of all physical data in a computational cell
class Variable {
 private:
  // disallow copy constructor
  Variable(Variable const &var) = delete;

 public:
  enum { Size = NHYDRO + NCLOUD + NTRACER + NCHEMISTRY + NSTATIC };
  enum class Type {
    Prim = 0,
    Cons = 1,
    Frac = 2,
    Chem = 3
  }

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
  Variable() : type_(Type::Prim) {
    {
      w = data_.data();
      c = w + NHYDRO;
      s = c + NCLOUD;
      q = s + NTRACER;
      x = q + NCHEMISTRY;
    }

    // Assignment operator
    Variable &operator=(const Variable &other) {
      // Check for self-assignment
      if (this == &other) {
        return *this;
      }

      // Perform member-wise assignment
      std::copy(other.w, other.w + Size, w);

      return *this;
    }

    void SetType(Type type) { type_ = type; }

    Type GetType() const { return type_; }

    void ConvertTo(Type type);

   private:
    // data holder
    std::array<Real, Size> data_;

    // type
    Type type_;
  };

#endif  // SRC_VARIABLE_HPP_
