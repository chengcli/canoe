#ifndef SRC_OPACITY_GIANTS_FREEDMAN_HPP_
#define SRC_OPACITY_GIANTS_FREEDMAN_HPP_

// C/C++
#include <string>
#include <vector>

// opacity
#include <opacity/absorber.hpp>

class ParameterInput;

// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  static const Real c1, c2, c3, c4, c5, c6, c7;

  FreedmanMean() : Absorber("FreedmanMean") {}
  virtual ~FreedmanMean() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class FreedmanMean2 : public Absorber {
 public:
  static const Real c1, c2, c3, c4, c5, c6, c7, c13;

  FreedmanMean2() : Absorber("FreedmanMean") {}
  virtual ~FreedmanMean2() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class FreedmanSimple : public Absorber {
 public:
  FreedmanSimple() : Absorber("FreedmanSimple") {}
  virtual ~FreedmanSimple() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

class FreedmanSimple2 : public Absorber {
 public:
  FreedmanSimple2() : Absorber("FreedmanSimple") {}
  virtual ~FreedmanSimple2() {}
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;
};

#endif  // SRC_OPACITY_GIANTS_FREEDMAN_HPP_
