#ifndef SRC_SNAP_MESHBLOCK_IMPL_HPP_
#define SRC_SNAP_MESHBLOCK_IMPL_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/mesh/mesh.hpp>

// snap
#include "decomposition/decomposition.hpp"
#include "implicit/implicit_solver.hpp"
#include "thermodynamics/thermodynamics.hpp"

class ParameterInput;
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

  ThermodynamicsPtr pthermo;
  DecompositionPtr pdec;
  ImplicitSolverPtr phevi;

  FaceReconstruct *precon;
  Forcing *pforce;

  Radiation *prad;
  std::vector<Inversion *> fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

#endif  // SRC_SNAP_MESHBLOCK_IMPL_HPP_
