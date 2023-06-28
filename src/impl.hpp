#ifndef SRC_IMPL_HPP_
#define SRC_IMPL_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/mesh/mesh.hpp>

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"

// harp
#include "harp/radiation.hpp"

// inversion
#include "inversion/inversion.hpp"

class ParameterInput;
class FaceReconstruct;
class Forcing;

//! \class MeshBlock::Impl
//  \brief opaque pointer class implements additional functionality of Athena
//  MeshBlock
class MeshBlock::Impl {
 public:
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl() {}

  AthenaArray<Real> du;  // stores tendency

  ThermodynamicsPtr pthermo;
  DecompositionPtr pdec;
  ImplicitSolverPtr phevi;

  // FaceReconstruct *precon;
  // Forcing *pforce;

  RadiationPtr prad;
  InversionQueue fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

#endif  // SRC_IMPL_HPP_
