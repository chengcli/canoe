#ifndef SRC_VARIABLE_HPP_
#define SRC_VARIABLE_HPP_

// C/C++
#include <array>

// canoe
#include <configure.hpp>

//! \class Variable
//  \brief a collection of all physical data in a computational cell
class Variable {
 private:
  // disallow copy constructor
  Variable(Variable const &var) = delete;

 public:
  enum { Size = NHYDRO + NCLOUD + NTRACER + NCHEMISTRY + NSTATIC };
  enum class Type {
    MassFrac = 0,
    MassConc = 1,
    MoleFrac = 2,
    MoleConc = 3
  }

  //! data pointers
  //! hydro data
  Real *w;

  //! cloud data
  Real *c;

  //! tracer data
  Real *x;

  //! chemistry data
  Real *q;

  //! static data
  Real *s;

  //! particle data
  Real *d;

  // constructor
  Variable() : type_(Type::Prim) {
    w = data_.data();
    c = w + NHYDRO;
    x = c + NCLOUD;
    q = x + NTRACER;
    s = q + NCHEMISTRY;
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

  void SetType(Type type) { mytype_ = type; }

  Type GetType() const { return mytype_; }

  void ConvertToPrimitive();
  void ConvertToConserved();
  void ConvertToMoleFraction();
  void ConvertToMoleConcentration();

 protected:
  void primitiveToConserved();
  void primitiveToMoleFraction();
  void primitiveToMoleConcentration();
  void conservedToPrimitive();
  void conservedToMoleFraction();
  void conservedToMoleConcentration();
  void moleFractionToPrimitive();
  void moleFractionToConserved();
  void moleFractionToMoleConcentration();
  void moleConcentrationToPrimitive();
  void moleConcentrationToConserved();
  void moleConcentrationToMoleFraction();

  // data holder
  std::array<Real, Size> data_;

  // type
  Type mytype_;
};

#endif  // SRC_VARIABLE_HPP_
