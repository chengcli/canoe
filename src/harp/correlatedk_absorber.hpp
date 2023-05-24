#ifndef SRC_HARP_CORRELATEDK_ABSORBER_HPP_
#define SRC_HARP_CORRELATEDK_ABSORBER_HPP_

// C/C++
#include <string>
#include <vector>

// harp
#include "absorber.hpp"

class CorrelatedKAbsorber : public Absorber {
 public:
  CorrelatedKAbsorber(MeshBlock *pmb, ParameterInput *pin, std::string bname,
                      std::string name)
      : Absorber(pmb, pin, bname, name) {}
  virtual ~CorrelatedKAbsorber() {}
  void loadCoefficient(std::string fname, int bid);
  // Real ckAbsorptionCoefficient(int mw, int mg, Real const prim[]) const;
  Real getAttenuation(Real g1, Real g2, CellVariables const &var) const;

 protected:
  int len_[3];               /**< length of interpolation axis */
  std::vector<Real> axis_;   /**< interpolation axis */
  std::vector<Real> kcoeff_; /**< absorption coefficient */
};

#endif  // SRC_HARP_CORRELATEDK_ABSORBER_HPP_
