#ifndef SRC_OPACITY_HITRAN_ABSORBER_HPP_
#define SRC_OPACITY_HITRAN_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include <harp/absorber.hpp>

class HitranAbsorber : public Absorber {
  friend std::ostream& operator<<(std::ostream& os, HitranAbsorber const& ab);
  // friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //   int const *ck_axis, Real const *ck_wave, int nbins);
 public:
  HitranAbsorber(std::string name, SpeciesNames const& species,
                 ParameterMap params)
      : Absorber(name, species, params) {}
  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname, size_t bid) override;
  Real GetAttenuation(Real wave1, Real wave2,
                      Variable const& var) const override;

 protected:
  size_t len_[3];            /**< length of interpolation axis */
  std::vector<Real> axis_;   /**< interpolation axis */
  std::vector<Real> kcoeff_; /**< absorption coefficient */
  AthenaArray<Real> refatm_; /**< reference atmosphere */
  Real getRefTemp(Real pres) const;
};

#endif  // SRC_HARP_HITRAN_ABSORBER_HPP_
