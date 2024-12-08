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
  Real GetAttenuation(int m, AirParcel const& var) const override;
  void ModifySpectralGrid(std::vector<SpectralBin>& spec) const override;

 protected:
  size_t len_[4];

  std::vector<Real> axis_;
  std::vector<Real> p_;
  std::vector<Real> t_;
  std::vector<Real> bin_centers_; 
  std::vector<Real> bin_edges_;
  std::vector<Real> kcoeff_;
  std::vector<Real> weights_;
};

class HeliosCK : public AbsorberCK {
 public:
  HeliosCK(std::string name);
  virtual ~HeliosCK() {}

  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real g1, Real g2, AirParcel const& var) const override;
  void ModifySpectralGrid(std::vector<SpectralBin>& spec) const override;
};

#endif  // SRC_OPACITY_ABSORBER_CK_HPP_
