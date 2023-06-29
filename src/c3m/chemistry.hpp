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
  Chemistry(MeshBlock *pmb, ParameterInput *pin);

 protected:
  AthenaArray<Real> r, s;
};

using ChemistryPtr = std::shared_ptr<Cloud>;

#endif  // SRC_C3M_CHEMISTRY_HPP_
