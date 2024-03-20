#ifndef SRC_C3M_CHEMISTRY_HPP_
#define SRC_C3M_CHEMISTRY_HPP_

// C/C++
#include <memory>

// athenapp
#include <athena/athena.hpp>

// cantera
#include <cantera/OneDim.h>

class MeshBlock;
class ParameterInput;

class AtmChemistry : public Cantera::OneDim {
 public:
  AthenaArray<Real> w, u;

  AtmChemistry(MeshBlock *pmb, ParameterInput *pin);
  ~AtmChemistry();

 protected:
  void fromMeshBlock(MeshBlock *pmb);
  MeshBlock *pmy_block_;
};

using AtmChemistryPtr = std::shared_ptr<AtmChemistry>;

#endif  // SRC_C3M_CHEMISTRY_HPP_
