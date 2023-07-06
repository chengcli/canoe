#ifndef SRC_C3M_CHEMISTRY_HPP_
#define SRC_C3M_CHEMISTRY_HPP_

// C/C++
#include <memory>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

class Chemistry {
 public:
  AthenaArray<Real> w, u;

  Chemistry(MeshBlock *pmb, ParameterInput *pin);
  ~Chemistry();

 protected:
  MeshBlock *pmy_block_;
};

using ChemistryPtr = std::shared_ptr<Chemistry>;

#endif  // SRC_C3M_CHEMISTRY_HPP_
