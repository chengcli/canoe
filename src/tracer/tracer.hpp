#ifndef SRC_TRACER_TRACER_HPP_
#define SRC_TRACER_TRACER_HPP_

// C/C++
#include <memory>
#include <vector>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;

class Tracer {
 public:
  AthenaArray<Real> r, s;

  Tracer(MeshBlock *pmb, ParameterInput *pin);

 protected:
};

using TracerPtr = std::shared_ptr<Tracer>;

#endif  //  SRC_TRACER_TRACER_HPP_
