#ifndef SRC_HARP_FREEDMAN_SIMPLE_HPP_
#define SRC_HARP_FREEDMAN_SIMPLE_HPP_

#include <string>

#include "absorber.hpp"

class FreedmanSimple : public Absorber {
 public:
  FreedmanSimple(MeshBlock *pmb, ParameterInput *pin, std::string bname);
  virtual ~FreedmanSimple() {}
  Real getAttenuation(Real wave1, Real wave2, CellVariables const &var) const;

 private:
  Real scale_;
};

#endif  // SRC_HARP_FREEDMAN_SIMPLE_HPP_
