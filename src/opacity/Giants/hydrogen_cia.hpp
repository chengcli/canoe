#ifndef SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP
#define SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP

// C/C++
#include <iostream>
#include <vector>

// harp
#include <harp/absorber.hpp>

class XizH2H2CIA : public Absorber {
 public:
  XizH2H2CIA(SpeciesNames const& species, ParameterMap params)
      : Absorber("H2-H2-CIA", species, params) {
    category_ = "cia";
  }

  virtual ~XizH2H2CIA() {}
  void LoadCoefficient(std::string fname, int bid = -1);
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class XizH2HeCIA : public Absorber {
 public:
  XizH2HeCIA(SpeciesNames const& species, ParameterMap params)
      : Absorber("H2-He-CIA", species, params) {
    category_ = "cia";
  }

  virtual ~XizH2HeCIA() {}
  void LoadCoefficient(std::string fname, int bid = -1);
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;

 protected:
  Real mixr2_;
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class OrtonCIA : public Absorber {
 public:
  OrtonCIA(SpeciesNames const& species, ParameterMap params)
      : Absorber("Orthon-CIA", species, params) {
    category_ = "cia";
  }

  virtual ~OrtonCIA() {}
  void LoadCoefficient(std::string fname, int bid = -1);
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

#endif
