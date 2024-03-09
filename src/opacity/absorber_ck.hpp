#ifndef SRC_OPACITY_ABSORBER_CK_HPP_
#define SRC_OPACITY_ABSORBER_CK_HPP_

// C/C++
#include <string>
#include <vector>

// opacity
#include "absorber.hpp"

class AbsorberCK : public Absorber {
 public:
  AbsorberCK(std::string name) : Absorber(name) {}
  virtual ~AbsorberCK() {}

  void LoadCoefficient(std::string fname, size_t bid) override;
  Real GetAttenuation(Real g1, Real g2, AirParcel const& var) const override;

 protected:
  //! shape of interpolation axes, nwave, npres, ntemp
  size_t len_[3];

  //! interpolation axes
  std::vector<Real> axis_;

  //! absorption coefficients
  std::vector<Real> kcoeff_;
};

class HeliosCKPremix : public AbsorberCK {
 public:
  HeliosCKPremix(std::string name) : AbsorberCK(name) {}
  virtual ~HeliosCKPremix() {}

  void LoadCoefficient(std::string fname, size_t bid) override;
  Real GetAttenuation(Real g1, Real g2, AirParcel const& var) const override;
};

#endif  // SRC_OPACITY_ABSORBER_CK_HPP_
