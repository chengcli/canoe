#ifndef SRC_FLASK_CHEMISTRY_HPP_
#define SRC_FLASK_CHEMISTRY_HPP_

// C/C++
#include <memory>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

class Flask {
 public:
  AthenaArray<Real> w, u;

  Flask(MeshBlock *pmb, ParameterInput *pin);
  ~Flask();

 protected:
  MeshBlock *pmy_block_;
};

using ChemistryPtr = std::shared_ptr<Flask>;

#endif  // SRC_C3M_CHEMISTRY_HPP_
