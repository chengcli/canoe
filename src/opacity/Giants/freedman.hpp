#ifndef SRC_OPACITY_GIANTS_FREEDMAN_HPP_
#define SRC_OPACITY_GIANTS_FREEDMAN_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include <harp/absorber.hpp>

class ParameterInput;

// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  FreedmanMean(std::vector<std::string> species, ParameterMap params)
      : Absorber("FreedmanMean", species, params) {}
  virtual ~FreedmanMean() {}
  Real GetAttenuation(Real wave1, Real wave2, Variable const &var) const;
};

class FreedmanSimple : public Absorber {
 public:
  FreedmanSimple(std::vector<std::string> species, ParameterMap params)
      : Absorber("FreedmanSimple", species, params) {}
  virtual ~FreedmanSimple() {}
  Real GetAttenuation(Real wave1, Real wave2, Variable const &var) const;
};

#endif  // SRC_OPACITY_GIANTS_FREEDMAN_HPP_
