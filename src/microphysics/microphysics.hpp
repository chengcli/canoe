#ifndef SRC_MICROPHYSICS_MICROPHYSICS_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICS_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <air_parcel.hpp>
#include <schedulers.hpp>
#include <virtual_groups.hpp>

class MeshBlock;
class ParameterInput;
class MicrophysicalSchemeBase;
class Hydro;

class Microphysics {
 public:
  // tem, v1, v2, v3
  // enum { NCLOUD_HYDRO = 4 };
  static const std::string input_key;

  // access members
  AthenaArray<Real> w, u, vsedf[3], mass_flux[3];

  /// constructor and destructor
  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();

  /// functions
  size_t GetNumSystems() const { return systems_.size(); }
  std::shared_ptr<MicrophysicalSchemeBase> GetSystem(int i) const {
    return systems_[i];
  }
  // void AddFrictionalHeating(std::vector<AirParcel> &air_column) const;
  void EvolveSystems(AirColumn &air_column, Real time, Real dt);

  /// inbound functions
  void SetSedimentationVelocityFromConserved(Hydro const *phydro);

 protected:
  AthenaArray<Real> vsed_[3];
  std::vector<std::shared_ptr<MicrophysicalSchemeBase>> systems_;

 private:
  MeshBlock const *pmy_block_;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;

namespace AllTasks {

// hydro tasks should be move into hydro in the future
bool hydro_implicit_correction(MeshBlock *pmb, IntegrationStage stage);
bool hydro_calculate_flux(MeshBlock *pmb, IntegrationStage stage);

bool microphysics_set_sedimentaton_velocity(MeshBlock *pmb,
                                            IntegrationStage stage);
bool microphysics_set_mass_flux(MeshBlock *pmb, IntegrationStage stage);
bool microphysics_evolve_system(MeshBlock *pmb, IntegrationStage stage);

// scalar tasks should be move into scalars in the future
bool scalar_calculate_flux(MeshBlock *pmb, IntegrationStage stage);

}  // namespace AllTasks

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
