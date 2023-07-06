#ifndef SRC_IMPL_HPP_
#define SRC_IMPL_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include "variable.hpp"

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"

// harp
#include "harp/radiation.hpp"

// cloud
#include "dusts/cloud.hpp"

// tracer
#include "tracer/tracer.hpp"

// tracer
#include "c3m/chemistry.hpp"

// inversion
#include "inversion/inversion.hpp"

class ParameterInput;

//! \class MeshBlock::Impl
//  \brief opaque pointer class implements additional functionality of Athena
//  MeshBlock
class MeshBlock::Impl {
 public:
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  AthenaArray<Real> du;  // stores tendency

  DecompositionPtr pdec;
  ImplicitSolverPtr phevi;

  CloudPtr pcloud;
  ChemistryPtr pchem;
  TracerPtr ptracer;
  // StaticVariablePtr pstatic;

  RadiationPtr prad;
  InversionQueue fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

  void GatherFromPrimitive(Variable *prim, int k, int j, int i);
  void GatherFromConserved(Variable *cons, int k, int j, int i);

  void DistributeToPrimitive(Variable const& prim, int k, int j, int i);
  void DistributeToConserved(Variable const& cons, int k, int j, int i);

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

int find_pressure_level_lesser(Real pmax, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie);

#endif  // SRC_IMPL_HPP_
