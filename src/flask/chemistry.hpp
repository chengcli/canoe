#ifndef SRC_FLASK_CHEMISTRY_HPP_
#define SRC_FLASK_CHEMISTRY_HPP_

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
};

using ChemistryPtr = std::shared_ptr<Chemistry>;

#endif  // SRC_FLASK_CHEMISTRY_HPP_
