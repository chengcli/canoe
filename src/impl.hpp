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

class MeshOutputGroup;
class FITSOutputGroup;

class Decomposition;
class ImplicitSolver;
class Tracer;
class TurbulenceModel;
class CelestrialBody;
class CubedSphere;
class SingleColumn;
// class Surface;

class Scheduler;

class Diagnostics;
class Forcing;

class ExchangerBase;

namespace harp {
class RadiationImpl;
}  // namespace harp

namespace kintera {
class KineticsImpl;
}  // namespace kintera

namespace snap {
class IdealMoistImpl;
class SedHydroImpl;
}  // namespace snap

//! \class MeshBlock::Impl
//! \brief opaque pointer class implements all physics on a MeshBlock
class MeshBlock::Impl {
 public:
  /// public data
  AthenaArray<Real> du;  // stores tendency

  std::shared_ptr<Decomposition> pdec;
  std::shared_ptr<ImplicitSolver> phevi;

  std::shared_ptr<snap::IdealMoistImpl> peos;
  std::shared_ptr<kintera::KineticsImpl> pkinet;
  std::shared_ptr<harp::RadiationImpl> prad;
  std::shared_ptr<snap::SedHydroImpl> psed;

  std::shared_ptr<TurbulenceModel> pturb;
  std::shared_ptr<CelestrialBody> planet;
  std::shared_ptr<CubedSphere> pexo3;
  std::shared_ptr<SingleColumn> pscm;

  std::vector<std::shared_ptr<Diagnostics>> all_diags;
  std::vector<std::shared_ptr<Forcing>> all_forcings;

 public:  // constructor and destructor
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

 public:  // member functions
  ExchangerBase *FindExchanger(char const *name) const {
    return exchangers_.at(name);
  }

  auto &GetMeshOutputGroups() const { return mesh_outputs_; }
  auto &GetFITSOutputGroups() const { return fits_outputs_; }

 protected:
  std::map<char const *, ExchangerBase *> exchangers_;

  std::vector<std::weak_ptr<MeshOutputGroup>> mesh_outputs_;
  std::vector<std::weak_ptr<FITSOutputGroup>> fits_outputs_;

 private:
  MeshBlock const *pmy_block_;
};

#endif  // SRC_IMPL_HPP_
