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
#include "clouds/cloud.hpp"

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

  TracerPtr ptracer;
  ChemistryPtr pchem;
  // StaticVariablePtr pstatic;

  CloudPtr pcloud;

  RadiationPtr prad;
  InversionQueue fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

  void GatherPrimitive(Variable *prim, int k, int j, int i);
  void GatherConserved(Variable *cons, int k, int j, int i);
  void GatherMoleFraction(Variable *qfrac, int k, int j, int i);
  void GatherMoleConcentration(Variable *qconc, int k, int j, int i);

  void DistributePrimitive(Variable &prim, int k, int j, int i);
  void DistributeConserved(Variable &cons, int k, int j, int i);

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

int find_pressure_level_lesser(Real pmax, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie);

#endif  // SRC_IMPL_HPP_
