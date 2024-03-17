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

  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real g1, Real g2, AirParcel const& var) const override;

 protected:
  //! shape of interpolation axes, ntemp, npres, ngpoints
  size_t len_[3];

  //! interpolation axes
  std::vector<Real> axis_;

  //! absorption coefficients
  std::vector<Real> kcoeff_;

  //! g-point weights
  std::vector<Real> weights_;
};

class HeliosCK : public AbsorberCK {
 public:
  HeliosCK(std::string name) : AbsorberCK(name) {}
  virtual ~HeliosCK() {}

  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real g1, Real g2, AirParcel const& var) const override;
  void ModifySpectralGrid(std::vector<SpectralBin>& spec) const override;
};

#endif  // SRC_OPACITY_ABSORBER_CK_HPP_
