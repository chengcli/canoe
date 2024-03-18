#ifndef SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP
#define SRC_OPACITY_GIANTS_HYDROGEN_CIA_HPP

// C/C++
#include <iostream>
#include <vector>

// opacity
#include <opacity/absorber.hpp>

class XizH2H2CIA : public Absorber {
 public:
  XizH2H2CIA();

  virtual ~XizH2H2CIA() {}
  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    Real k1 = getAttenuation1(wave1, var);
    Real k2 = getAttenuation1(wave2, var);
    return (k1 + k2) / 2.;
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class XizH2HeCIA : public Absorber {
 public:
  XizH2HeCIA();

  virtual ~XizH2HeCIA() {}
  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    Real k1 = getAttenuation1(wave1, var);
    Real k2 = getAttenuation1(wave2, var);
    return (k1 + k2) / 2.;
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;

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
  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

#endif
