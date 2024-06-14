#ifndef SRC_OPACITY_HITRAN_ABSORBER_HPP_
#define SRC_OPACITY_HITRAN_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// opacity
#include "absorber.hpp"

class HitranAbsorber : public Absorber {
  friend std::ostream& operator<<(std::ostream& os, HitranAbsorber const& ab);
  // friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //   int const *ck_axis, Real const *ck_wave, int nbins);
 public:
  HitranAbsorber(std::string name) : Absorber(name) {}

  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname, int bid) override;
  Real GetAttenuation(Real wave1, Real wave2,
                      AirParcel const& var) const override;

 protected:
  size_t len_[3];            /**< length of interpolation axis */
  std::vector<Real> axis_;   /**< interpolation axis */
  std::vector<Real> kcoeff_; /**< absorption coefficient */
  AthenaArray<Real> refatm_; /**< reference atmosphere */
  Real getRefTemp(Real pres) const;
};

class HitranAbsorberCK : public HitranAbsorber {
 public:
  using Base = HitranAbsorber;
  HitranAbsorberCK(std::string name) : HitranAbsorber(name) {}

  void LoadCoefficient(std::string fname, int bid) override;
  void ModifySpectralGrid(std::vector<SpectralBin>& spec) const override;

 private:
  std::vector<Real> weights_;
};

#endif  // SRC_HARP_HITRAN_ABSORBER_HPP_
