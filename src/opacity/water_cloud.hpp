#ifndef SRC_OPACITY_WATER_CLOUD_HPP_
#define SRC_OPACITY_WATER_CLOUD_HPP_

// C/C++
#include <string>
#include <vector>

// canoe
#include <air_parcel.hpp>

// harp
#include "absorber.hpp"

class FuWaterLiquidCloud : public Absorber {
 public:
  FuWaterLiquidCloud() : Absorber("H2O(l)") {}

  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2, var);
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 AirParcel const& var) const override {
    return getSingleScatteringAlbedo1((wave1 + wave2) / 2., var);
  }

  void GetPhaseMomentum(Real* pp, Real wave1, Real wave2, AirParcel const& var,
                        int np) const override {
    getPhaseMomentum1(pp, (wave1 + wave2) / 2., var, np);
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;
  Real getSingleScatteringAlbedo1(Real wave, AirParcel const& var) const;
  void getPhaseMomentum1(Real* pp, Real wave, AirParcel const& var,
                         int np) const;
};

class FuWaterIceCloud : public Absorber {
 public:
  FuWaterIceCloud() : Absorber("H2O(s)") {}

  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2, var);
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 AirParcel const& var) const override {
    return getSingleScatteringAlbedo1((wave1 + wave2) / 2., var);
  }

  void GetPhaseMomentum(Real* pp, Real wave1, Real wave2, AirParcel const& var,
                        int np) const override {
    getPhaseMomentum1(pp, (wave1 + wave2) / 2., var, np);
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;
  Real getSingleScatteringAlbedo1(Real wave, AirParcel const& var) const;
  void getPhaseMomentum1(Real* pp, Real wave, AirParcel const& var,
                         int np) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

class XuWaterIceCloud : public Absorber {
 public:
  XuWaterIceCloud() : Absorber("H2O(s)") {}

  void LoadCoefficient(std::string fname, int bid) override;

  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2, var);
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 AirParcel const& var) const override {
    return getSingleScatteringAlbedo1((wave1 + wave2) / 2., var);
  }

  void GetPhaseMomentum(Real* pp, Real wave1, Real wave2, AirParcel const& var,
                        int np) const override {
    getPhaseMomentum1(pp, (wave1 + wave2) / 2., var, np);
  }

 protected:
  Real getAttenuation1(Real wave, AirParcel const& var) const;
  Real getSingleScatteringAlbedo1(Real wave, AirParcel const& var) const;
  void getPhaseMomentum1(Real* pp, Real wave, AirParcel const& var,
                         int np) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

#endif  // SRC_OPACITY_WATER_CLOUD_HPP_
