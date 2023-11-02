#ifndef SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP
#define SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP

// C/C++
#include <iostream>
#include <vector>

// opacity
#include <opacity/absorber.hpp>

class XizH2H2CIA : public Absorber {
 public:
  XizH2H2CIA() : Absorber("H2-H2-CIA") {}

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
  XizH2HeCIA() : Absorber("H2-He-CIA") {}

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
  OrtonCIA() : Absorber("Orthon-CIA") {}

  virtual ~OrtonCIA() {}
  void LoadCoefficient(std::string fname, int bid = -1);
  Real GetAttenuation(Real wave1, Real wave2, AirParcel const& var) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

#endif
