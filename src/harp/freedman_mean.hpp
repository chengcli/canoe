#ifndef SRC_HARP_FREEDMAN_MEAN_HPP_
#define SRC_HARP_FREEDMAN_MEAN_HPP_

// Athena++ header
#include <string>
#include "absorber.hpp"

// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  FreedmanMean(MeshBlock *pmb, ParameterInput *pin, std::string bname)
      : Absorber(pmb, pin, bname, "FreedmanMean") {}
  virtual ~FreedmanMean() {}
  Real getAttenuation(Real wave1, Real wave2, CellVariables const &var) const;
};

#endif  // SRC_HARP_FREEDMAN_MEAN_HPP_
