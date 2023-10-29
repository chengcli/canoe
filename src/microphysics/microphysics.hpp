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

// microphysics
#include "microphysical_scheme.hpp"

class MeshBlock;
class ParameterInput;

class Microphysics {
 public:
  // tem, v1, v2, v3
  // enum { NCLOUD_HYDRO = 4 };

  // access members
  AthenaArray<Real> w, u;

  /// constructor and destructor
  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();

  /// functions
  size_t GetNumSystems() const { return systems_.size(); }
  std::string GetSystemName(int i) const { return systems_[i]->GetName(); }
  std::shared_ptr<MicrophysicalScheme> GetSystem(int i) const {
    return systems_[i];
  }

  // void AddFrictionalHeating(std::vector<AirParcel> &air_column) const;

  void EvolveSystems(AirColumn &air_column, Real time, Real dt);

  AthenaArray<Real> const &GetSedimentationVelocityFace(int dir) const {
    return vsedf_[dir];
  }

  // set by Riemann Solver
  void SetFluidMassFlux(int dir, int n, int k, int j, int i, Real val) {
    mass_flux_[dir](n, k, j, i) = val;
  }

  AthenaArray<Real> const &GetFluidMassFlux(int dir) const {
    return mass_flux_[dir];
  }

  // called in tasklist UpdateAllConserved
  void UpdateSedimentationVelocityFromConserved();

 protected:
  AthenaArray<Real> mass_flux_[3];
  AthenaArray<Real> vsed_[3], vsedf_[3];
  std::vector<MicrophysicalSchemePtr> systems_;

 private:
  MeshBlock const *pmy_block_;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
