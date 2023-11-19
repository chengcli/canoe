#ifndef SRC_OPACITY_SIMPLE_CLOUD_HPP_
#define SRC_OPACITY_SIMPLE_CLOUD_HPP_

// C/C++
#include <string>
#include <vector>

// canoe
#include <air_parcel.hpp>

// harp
#include "absorber.hpp"

class SimpleCloud : public Absorber {
 public:
  SimpleCloud(std::string name) : Absorber(name) {}

  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override {
    Real a = getAttenuation1(wave1, var);
    Real b = getAttenuation1(wave2, var);
    return (a + b) / 2.;
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 AirParcel const& var) const override {
    Real a = getSingleScatteringAlbedo1(wave1, var);
    Real b = getSingleScatteringAlbedo1(wave2, var);
    return (a + b) / 2.;
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

#endif  // SRC_OPACITY_SIMPLE_CLOUD_HPP_
