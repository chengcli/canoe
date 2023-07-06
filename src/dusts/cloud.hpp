#ifndef SRC_DUSTS_CLOUD_HPP_
#define SRC_DUSTS_CLOUD_HPP_

// C/C++
#include <memory>
#include <vector>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

class Cloud {
 public:
  AthenaArray<Real> w, u;

  Cloud(MeshBlock *pmb, ParameterInput *pin);
  ~Cloud();

 protected:
  MeshBlock *pmy_block_;
};

using CloudPtr = std::shared_ptr<Cloud>;

#endif  // SRC_CLOUDS_CLOUD_HPP_
