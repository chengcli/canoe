#ifndef SRC_VARIABLE_HPP_
#define SRC_VARIABLE_HPP_

// C/C++
#include <array>
#include <iostream>

// athena
#include <athena/athena.hpp>  // Real

// canoe
#include <configure.hpp>

//! \class Variable
//  \brief a collection of all physical data in a computational cell
class Variable {
 public:
  enum { Size = NHYDRO + NCLOUD + NCHEMISTRY + NTRACER + NSTATIC };

  enum class Type { MassFrac = 0, MassConc = 1, MoleFrac = 2, MoleConc = 3 };

  friend std::ostream &operator<<(std::ostream &os, Type const &type);
  friend std::ostream &operator<<(std::ostream &os, Variable const &var);

  //! data pointers
  //! hydro data
  Real *w;

  //! cloud data
  Real *c;

  //! chemistry data
  Real *q;

  //! tracer data
  Real *x;

  //! static data
  Real *s;

  //! particle data
  Real *d;

  // constructor
  explicit Variable(Type type = Type::MoleFrac) : mytype_(type) {
    w = data_.data();
    c = w + NHYDRO;
    q = c + NCLOUD;
    x = q + NCHEMISTRY;
    s = x + NTRACER;
    d = s + NSTATIC;
    std::fill(data_.begin(), data_.end(), 0.0);
  }

  // copy constructor
  Variable(Variable const &other) : mytype_(other.mytype_) {
    data_ = other.data_;
    w = data_.data();
    c = w + NHYDRO;
    q = c + NCLOUD;
    x = q + NCHEMISTRY;
    s = x + NTRACER;
    d = s + NSTATIC;
  }

  // Assignment operator
  Variable &operator=(const Variable &other) {
    // Check for self-assignment
    if (this == &other) {
      return *this;
    }

    data_ = other.data_;
    return *this;
  }

  void SetType(Type type) { 
    mytype_ = type; 
  }

  Type GetType() const { return mytype_; }

  void SetZero() { std::fill(data_.begin(), data_.end(), 0.0); }

  void ConvertTo(Variable::Type type);

  void ConvertToMassFraction();
  void ConvertToMassConcentration();
  void ConvertToMoleFraction();
  void ConvertToMoleConcentration();

 protected:
  void massFractionToMassConcentration();
  void massConcentrationToMassFraction();

  void massFractionToMoleFraction();
  void moleFractionToMassFraction();

  void massConcentrationToMoleFraction();
  void moleFractionToMassConcentration();

  void moleFractionToMoleConcentration();
  void moleConcentrationToMoleFraction();

  void massConcentrationToMoleConcentration();
  void moleConcentrationToMassConcentration();

  void massFractionToMoleConcentration();
  void moleConcentrationToMassFraction();

  // data holder
  std::array<Real, Size> data_;

  // type
  Type mytype_;
};

#endif  // SRC_VARIABLE_HPP_
