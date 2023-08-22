#ifndef SRC_MICROPHYSICS_MICROPHYSICS_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICS_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// canoe
#include <air_parcel.hpp>

// microphysics
#include "microphysical_scheme.hpp"

class MeshBlock;
class ParameterInput;

class Microphysics {
 public:
  // tem, v1, v2, v3
  enum { NCLOUD_HYDRO = 4 };

  // access members
  AthenaArray<Real> w, u;

  Microphysics(MeshBlock *pmb, ParameterInput *pin);

  ~Microphysics();

  size_t GetNumSystems() const { return systems_.size(); }

  std::string GetSystemName(int i) const { return systems_[i]->GetName(); }

  MicrophysicalSchemePtr GetSystem(int i) const { return systems_[i]; }

  void AddFrictionalHeating(std::vector<AirParcel> &air_column) const;

  void EvolveSystems(std::vector<AirParcel> &air_column, Real time, Real dt);

  void SetSedimentationVelocity(int k, int j, int il, int iu);

  void AddSedimentationFlux(AthenaArray<Real> &sflx, int k, int j, int il,
                            int iu) const;

 protected:
  AthenaArray<Real> vsed_, vsedf_;

  AthenaArray<Real> hydro_;

  std::vector<MicrophysicalSchemePtr> systems_;

  MeshBlock *pmy_block_ = nullptr;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
