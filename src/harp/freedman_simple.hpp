#ifndef SRC_HARP_FREEDMAN_SIMPLE_HPP_
#define SRC_HARP_FREEDMAN_SIMPLE_HPP_

#include <string>

#include "absorber.hpp"

class MeshBlock;

class FreedmanSimple : public Absorber {
 public:
  FreedmanSimple(MeshBlock *pmb, ParameterInput *pin, std::string bname);
  virtual ~FreedmanSimple() {}
  Real GetAttenuation(Real wave1, Real wave2, CellVariables const &var) const;

 protected:
  MeshBlock *pmy_block_;
  Real scale_;
};

#endif  // SRC_HARP_FREEDMAN_SIMPLE_HPP_
