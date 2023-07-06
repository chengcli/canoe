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

 protected:
  Cloud(MeshBlock *pmb, ParameterInput *pin);
};

using CloudPtr = std::unique_ptr<Cloud>;
using CloudQueue = std::vector<CloudPtr>;

#endif  // SRC_CLOUDS_CLOUD_HPP_
