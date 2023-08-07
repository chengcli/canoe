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

  void GatherFromPrimitive(Variable *prim, int k, int j, int i) const;
  void GatherFromPrimitive(std::vector<Variable> &air_column, int k, int j,
                           int il, int iu) const {
    air_column.resize(iu - il + 1);
    for (int i = il; i <= iu; ++i) {
      GatherFromPrimitive(&air_column[i - il], k, j, i);
    }
  }

  void GatherFromConserved(Variable *cons, int k, int j, int i) const;
  void GatherFromConserved(std::vector<Variable> &air_column, int k, int j,
                           int il, int iu) const {
    air_column.resize(iu - il + 1);
    for (int i = il; i <= iu; ++i) {
      GatherFromConserved(&air_column[i - il], k, j, i);
    }
  }

  void DistributeToPrimitive(Variable const &prim, int k, int j, int i);
  void DistributeToPrimitive(std::vector<Variable> const &air_column, int k,
                             int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      DistributeToPrimitive(air_column[i - il], k, j, i);
    }
  }

  void DistributeToConserved(Variable const &cons, int k, int j, int i);
  void DistributeToConserved(std::vector<Variable> const &air_column, int k,
                             int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      DistributeToConserved(air_column[i - il], k, j, i);
    }
  }

  // TODO(cli) : more needs to be changed
  // called in task_list/time_integration.cpp
  void MapScalarsConserved(AthenaArray<Real> &s) {
    if (NCLOUD > 0) pcloud->u.InitWithShallowSlice(s, 4, 0, NCLOUD);
  }

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

int find_pressure_level_lesser(Real pmax, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie);

#endif  // SRC_IMPL_HPP_
