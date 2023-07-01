#ifndef SRC_CLOUDS_CLOUD_HPP_
#define SRC_CLOUDS_CLOUD_HPP_

// C/C++
#include <memory>
#include <vector>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

class Cloud {
 public:
  static std::vector<std::unique_ptr<Cloud>> NewCloudsQueue(
      MeshBlock *pmb, ParameterInput *pin);

 protected:
  Cloud(MeshBlock *pmb, ParameterInput *pin);
  AthenaArray<Real> r, s;
};

using CloudPtr = std::unique_ptr<Cloud>;
using CloudQueue = std::vector<CloudPtr>;

#endif  // SRC_CLOUDS_CLOUD_HPP_
