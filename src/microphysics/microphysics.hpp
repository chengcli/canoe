#ifndef SRC_MICROPHYSICS_MICROPHYSICS_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICS_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// snap
#include <air_parcel.hpp>

// microphysics
#include "microphysical_scheme.hpp"

class MeshBlock;
class ParameterInput;

class Microphysics {
 public:
  // tem, v1, v2, v3
  enum { NCLOUD_HYDRO = 4 };

  AthenaArray<Real> w, u;

  Microphysics(MeshBlock *pmb, ParameterInput *pin);

  ~Microphysics();

  void EvolveSystems(std::vector<AirParcel> &air_column, Real time, Real dt) {
    for (auto &system : systems_)
      for (auto &air : air_column) {
        system->AssembleReactionMatrix(system->GetRatePtr(),
                                       system->GetJacobianPtr(), air, time);
        system->EvolveOneStep(&air, time, dt);
      }
  }

  void AddFrictionalHeating(std::vector<AirParcel> &air_column);

  void AddSedimentationFlux(AthenaArray<Real> &u, int k, int j, int il, int iu);

 protected:
  AthenaArray<Real> vsed_, vsedf_;

  AthenaArray<Real> hydro_;

  std::vector<MicrophysicalSchemePtr> systems_;

  MeshBlock *pmy_block_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
