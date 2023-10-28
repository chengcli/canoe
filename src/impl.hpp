#ifndef SRC_IMPL_HPP_
#define SRC_IMPL_HPP_

// C/C++
#include <map>
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

class ParameterInput;

class Inversion;
class ParticleBase;

class MeshOutputGroup;
class FITSOutputGroup;

//! \class MeshBlock::Impl
//  \brief opaque pointer class implements additional functionality of Athena
//  MeshBlock
class MeshBlock::Impl {
 public:
  /// public data
  AthenaArray<Real> du;  // stores tendency

  DecompositionPtr pdec;
  ImplicitSolverPtr phevi;

  MicrophysicsPtr pmicro;
  ChemistryPtr pchem;
  TracerPtr ptracer;
  // StaticVariablePtr pstatic;

  TurbulenceModelPtr pturb;

  RadiationPtr prad;

  // std::vector<std::shared_ptr<Inversion>> all_fits;
  std::vector<std::shared_ptr<ParticleBase>> all_particles;

  /// constructor and destructor
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  /// functions
  void *GetExchanger(std::string name) const {
    return boundary_exchanger_.at(name);
  }

  auto GetMeshOutputGroups() const { return mesh_outputs_; }
  auto GetFITSOutputGroups() const { return fits_outputs_; }

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

  // TODO(cli) : more needs to be changed
  // called in task_list/time_integration.cpp
  void MapScalarsConserved(AthenaArray<Real> &s) {
    if (NCLOUD > 0) pmicro->u.InitWithShallowSlice(s, 4, 0, NCLOUD);
  }

 protected:
  std::map<std::string, void *> boundary_exchanger_;

  std::vector<std::weak_ptr<MeshOutputGroup>> mesh_outputs_;
  std::vector<std::weak_ptr<FITSOutputGroup>> fits_outputs_;

  Real reference_pressure_;
  Real pressure_scale_height_;

 private:
  MeshBlock *pmy_block_;
  int my_stage_;
};

// helper functions
int find_pressure_level_lesser(Real pmax, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie);

#endif  // SRC_IMPL_HPP_
