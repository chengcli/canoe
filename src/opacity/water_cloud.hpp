#ifndef SRC_OPACITY_WATER_CLOUD_HPP_
#define SRC_OPACITY_WATER_CLOUD_HPP_

// C/C++
#include <string>
#include <vector>

// canoe
#include <variable.hpp>

// harp
#include <harp/absorber.hpp>

class FuWaterLiquidCloud : public Absorber {
 public:
  FuWaterLiquidCloud(SpeciesNames const& species, ParameterMap params)
      : Absorber("H2O(l)", species, params) {}

  Real GetAttenuation(Real wave1, Real wave2,
                      Variable const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2, var);
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 Variable const& var) const override {
    return getSingleScatteringAlbedo1((wave1 + wave2) / 2., var);
  }

  void GetPhaseMomentum(Real* pp, Real wave1, Real wave2, Variable const& var,
                        int np) const override {
    getPhaseMomentum1(pp, (wave1 + wave2) / 2., var, np);
  }

 protected:
  Real getAttenuation1(Real wave, Variable const& var) const;
  Real getSingleScatteringAlbedo1(Real wave, Variable const& var) const;
  void getPhaseMomentum1(Real* pp, Real wave, Variable const& var,
                         int np) const;
};

class FuWaterIceCloud : public Absorber {
 public:
  FuWaterIceCloud(SpeciesNames const& species, ParameterMap params)
      : Absorber("H2O(s)", species, params) {}

  Real GetAttenuation(Real wave1, Real wave2,
                      Variable const& var) const override {
    return getAttenuation1((wave1 + wave2) / 2, var);
  }

  Real GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                 Variable const& var) const override {
    return getSingleScatteringAlbedo1((wave1 + wave2) / 2., var);
  }

  void GetPhaseMomentum(Real* pp, Real wave1, Real wave2, Variable const& var,
                        int np) const override {
    getPhaseMomentum1(pp, (wave1 + wave2) / 2., var, np);
  }

 protected:
  Real getAttenuation1(Real wave, Variable const& var) const;
  Real getSingleScatteringAlbedo1(Real wave, Variable const& var) const;
  void getPhaseMomentum1(Real* pp, Real wave, Variable const& var,
                         int np) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

class XuWaterIceCloud : public Absorber {
 public:
  XuWaterIceCloud(SpeciesNames const& species, ParameterMap params)
      : Absorber("H2O(s)", species, params) {}

  void LoadCoefficient(std::string fname, size_t bid) override;

 protected:
  Real getAttenuation1(Real wave, Variable const& var) const;
  Real getSingleScateringAlbedo1(Real wave, Variable const& var) const;
  void getPhaseMomentum(Real* pp, Real wave, Variable const& var, int np) const;

 protected:
  size_t len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

#endif  // SRC_OPACITY_WATER_CLOUD_HPP_
