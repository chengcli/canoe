#ifndef FREEDMAN_MEAN_HPP
#define FREEDMAN_MEAN_HPP

// Athena++ header
#include "absorber.hpp"

// Richard S. Freedman 2011. APJS
class FreedmanMean : public Absorber {
 public:
  FreedmanMean(MeshBlock *pmb, ParameterInput *pin, std::string bname)
      : Absorber(pmb, pin, bname, "FreedmanMean") {}
  virtual ~FreedmanMean() {}
  Real getAttenuation(Real wave1, Real wave2, CellVariables const &var) const;
};

#endif
