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

class ParameterInput;

class Decomposition;
class ImplicitSolver;
class Microphysics;
class Radiation;
class Chemistry;
class Tracer;
class TurbulenceModel;

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

  std::shared_ptr<Decomposition> pdec;
  std::shared_ptr<ImplicitSolver> phevi;
  std::shared_ptr<Microphysics> pmicro;
  std::shared_ptr<Radiation> prad;
  std::shared_ptr<Chemistry> pchem;
  std::shared_ptr<Tracer> ptracer;
  std::shared_ptr<TurbulenceModel> pturb;

  // StaticVariablePtr pstatic;

  std::vector<std::shared_ptr<Inversion>> all_fits;
  std::vector<std::shared_ptr<ParticleBase>> all_particles;

  /// constructor and destructor
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  /// functions
  void *GetExchanger(std::string name) const {
    return boundary_exchanger_.at(name);
  }

  auto &GetMeshOutputGroups() const { return mesh_outputs_; }
  auto &GetFITSOutputGroups() const { return fits_outputs_; }

  // TODO(cli) : more needs to be changed
  // called in task_list/time_integration.cpp
  void MapScalarsConserved(AthenaArray<Real> &s);

 protected:
  std::map<std::string, void *> boundary_exchanger_;
  std::vector<std::weak_ptr<MeshOutputGroup>> mesh_outputs_;
  std::vector<std::weak_ptr<FITSOutputGroup>> fits_outputs_;

 private:
  MeshBlock const *pmy_block_;
  int my_stage_;
};

#endif  // SRC_IMPL_HPP_
