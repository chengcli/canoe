#ifndef MESHBLOCK_IMPL_HPP
#define MESHBLOCK_IMPL_HPP

// C/C++ header
#include <memory>
#include <vector>

// Athena++ header
#include <athena/mesh/mesh.hpp>

class ParameterInput;
class Thermodynamics;
class Decomposition;
class ImplicitSolver;
class FaceReconstruct;
class Forcing;
class Radiation;
class Inversion;

//! \class MeshBlock::Impl
//  \brief opaque pointer class implements additional functionality of Athena
//  MeshBlock
class MeshBlock::Impl {
 public:
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  AthenaArray<Real> du;  // stores tendency

  Thermodynamics *pthermo;
  Decomposition *pdec;
  ImplicitSolver *phevi;
  FaceReconstruct *precon;
  Forcing *pforce;

  Radiation *prad;
  std::vector<Inversion *> fitq;

 private:
  MeshBlock const *pmy_block_;
};

#endif
