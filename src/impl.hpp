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
#include "air_parcel.hpp"

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"

// harp
#include "harp/radiation.hpp"

// microphysics
#include "microphysics/microphysics.hpp"

// tracer
#include "tracer/tracer.hpp"

// chemistry
#include "c3m/chemistry.hpp"

// turbulence
#include "snap/turbulence/turbulence_model.hpp"

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

  MicrophysicsPtr pmicro;
  ChemistryPtr pchem;
  TracerPtr ptracer;
  // StaticVariablePtr pstatic;

  TurbulenceModelPtr pturb;

  RadiationPtr prad;
  InversionQueue fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

  void GatherFromPrimitive(AirParcel *prim, int k, int j, int i) const;
  void GatherFromPrimitive(std::vector<AirParcel> &air_column, int k, int j,
                           int il, int iu) const {
    air_column.resize(iu - il + 1);
    for (int i = il; i <= iu; ++i) {
      GatherFromPrimitive(&air_column[i - il], k, j, i);
    }
  }

  void GatherFromConserved(AirParcel *cons, int k, int j, int i) const;
  void GatherFromConserved(std::vector<AirParcel> &air_column, int k, int j,
                           int il, int iu) const {
    air_column.resize(iu - il + 1);
    for (int i = il; i <= iu; ++i) {
      GatherFromConserved(&air_column[i - il], k, j, i);
    }
  }

  void DistributeToPrimitive(AirParcel const &prim, int k, int j, int i);
  void DistributeToPrimitive(std::vector<AirParcel> const &air_column, int k,
                             int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      DistributeToPrimitive(air_column[i - il], k, j, i);
    }
  }

  void DistributeToConserved(AirParcel const &cons, int k, int j, int i);
  void DistributeToConserved(std::vector<AirParcel> const &air_column, int k,
                             int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      DistributeToConserved(air_column[i - il], k, j, i);
    }
  }

  // TODO(cli) : more needs to be changed
  // called in task_list/time_integration.cpp
  void MapScalarsConserved(AthenaArray<Real> &s) {
    if (NCLOUD > 0) pmicro->u.InitWithShallowSlice(s, 4, 0, NCLOUD);
  }

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

int find_pressure_level_lesser(Real pmax, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie);

#endif  // SRC_IMPL_HPP_
