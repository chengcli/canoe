#ifndef SRC_HARP_FREEDMAN_MEAN_HPP_
#define SRC_HARP_FREEDMAN_MEAN_HPP_

// C/C++
#include <string>

// harp
#include "absorber.hpp"

class MeshBlock;
// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  FreedmanMean(MeshBlock *pmb, ParameterInput *pin, std::string bname)
      : Absorber("FreedmanMean") {}
  virtual ~FreedmanMean() {}
  Real GetAttenuation(Real wave1, Real wave2, Variable const &var) const;

 protected:
  MeshBlock *pmy_block_;
};

#endif  // SRC_HARP_FREEDMAN_MEAN_HPP_
