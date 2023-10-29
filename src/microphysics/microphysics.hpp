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

class MeshBlock;
class ParameterInput;
class MicrophysicalSchemeBase;

class Microphysics {
 public:
  // tem, v1, v2, v3
  // enum { NCLOUD_HYDRO = 4 };
  static const std::string input_key;

  // access members
  AthenaArray<Real> w, u, vsedf[3];

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

  /// outbound functions
  // 1. used by Riemann Solver, e.g. lmars.cpp
  void SetFluidMassFlux(int dir, int n, int k, int j, int i, Real val) {
    mass_flux_[dir](n, k, j, i) = val;
  }
  // 2. used by tracer, e.g. calculate_scalar_fluxes.cpp
  AthenaArray<Real> const &GetFluidMassFlux(int dir) const {
    return mass_flux_[dir];
  }
  // 3. used by tasklist UpdateAllConserved
  void UpdateSedimentationVelocityFromConserved();

 protected:
  AthenaArray<Real> mass_flux_[3];
  AthenaArray<Real> vsed_[3];
  std::vector<std::shared_ptr<MicrophysicalSchemeBase>> systems_;

 private:
  MeshBlock const *pmy_block_;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
