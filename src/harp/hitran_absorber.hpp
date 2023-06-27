#ifndef SRC_HARP_HITRAN_ABSORBER_HPP_
#define SRC_HARP_HITRAN_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include "absorber.hpp"

class HitranAbsorber : public Absorber {
  friend std::ostream& operator<<(std::ostream& os, HitranAbsorber const& ab);
  // friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //   int const *ck_axis, Real const *ck_wave, int nbins);
 public:
  HitranAbsorber(MeshBlock* pmb, ParameterInput* pin, std::string bname,
                 std::string name, int imol);
  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname, int bid);
  Real GetAttenuation(Real wave1, Real wave2, CellVariables const& var) const;

 protected:
  size_t len_[3];            /**< length of interpolation axis */
  std::vector<Real> axis_;   /**< interpolation axis */
  std::vector<Real> kcoeff_; /**< absorption coefficient */
  AthenaArray<Real> refatm_; /**< reference atmosphere */
  Real RefTemp_(Real pres) const;
};

#endif  // SRC_HARP_HITRAN_ABSORBER_HPP_
