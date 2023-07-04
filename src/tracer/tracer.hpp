#ifndef SRC_TRACER_TRACER_HPP_
#define SRC_TRACER_TRACER_HPP_

// C/C++
#include <memory>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

class Tracer {
 public:
  AthenaArray<Real> w, u;

  Tracer(MeshBlock *pmb, ParameterInput *pin);

 protected:
  MeshBlock *pmy_block_;
};

using TracerPtr = std::shared_ptr<Tracer>;

#endif  //  SRC_TRACER_TRACER_HPP_
