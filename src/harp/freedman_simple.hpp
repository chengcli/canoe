#ifndef FREEDMAN_SIMPLE_HPP
#define FREEDMAN_SIMPLE_HPP

// Athena++ header
#include "absorber.hpp"

class FreedmanSimple : public Absorber {
 public:
  FreedmanSimple(MeshBlock *pmb, ParameterInput *pin, std::string bname);
  virtual ~FreedmanSimple() {}
  Real getAttenuation(Real wave1, Real wave2, CellVariables const &var) const;

 private:
  Real scale_;
};

#endif
