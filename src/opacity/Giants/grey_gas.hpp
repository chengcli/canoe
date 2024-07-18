#ifndef SRC_OPACITY_GIANTS_GREY_GAS_HPP_
#define SRC_OPACITY_GIANTS_GREY_GAS_HPP_

// C/C++
#include <string>
#include <vector>

// opacity
#include <opacity/absorber.hpp>

class ParameterInput;

// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  static const Real c1, c2, c3, c4, c5, c6, c7, c13;

  FreedmanMean(std::string name) : Absorber(name) {}
  virtual ~FreedmanMean() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

//user defined
class SimpleGrey : public Absorber {
 public:
  SimpleGrey(std::string name) : Absorber(name) {}
  virtual ~SimpleGrey() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class JupGasv : public Absorber {
 public:
  JupGasv(std::string name) : Absorber(name) {}
  virtual ~JupGasv() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class JupGasir : public Absorber {
 public:
  JupGasir(std::string name) : Absorber(name) {}
  virtual ~JupGasir() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

#endif  // SRC_OPACITY_GIANTS_GREY_GAS_HPP_
