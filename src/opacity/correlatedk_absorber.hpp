#ifndef SRC_OPACITY_CORRELATEDK_ABSORBER_HPP_
#define SRC_OPACITY_CORRELATEDK_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include <harp/absorber.hpp>

class CorrelatedKAbsorber : public Absorber {
 public:
  CorrelatedKAbsorber(std::string name, SpeciesNames const& species,
                      ParameterMap params)
      : Absorber(name, species, params) {}
  virtual ~CorrelatedKAbsorber() {}
  void LoadCoefficient(std::string fname, size_t bid) override;
  // Real ckAbsorptionCoefficient(int mw, int mg, Real const prim[]) const;
  Real GetAttenuation(Real g1, Real g2, Variable const& var) const override;

 protected:
  size_t len_[3];            /**< length of interpolation axis */
  std::vector<Real> axis_;   /**< interpolation axis */
  std::vector<Real> kcoeff_; /**< absorption coefficient */
};

#endif  // SRC_OPACITY_CORRELATEDK_ABSORBER_HPP_
