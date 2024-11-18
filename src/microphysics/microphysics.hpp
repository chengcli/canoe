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
#include <virtual_groups.hpp>

// microphysics
#include "sedimentation.hpp"

class MeshBlock;
class ParameterInput;
class MicrophysicalSchemeBase;
class Hydro;

//! \brief root-level management class for microphysics
class Microphysics {
 public:  /// public access members
  //! \todo(CLI) track cloud temperature and momentum
  //! tem, v1, v2, v3

  //! microphysics input key in the input file [microphysics_config]
  static const std::string input_key;

  //! sedimentation velocity at cell interface [m/s]
  AthenaArray<Real> vsedf[3];

  //! mass flux of the dry fluid [kg/m^2/s]
  AthenaArray<Real> mass_flux[3];

 public:  /// constructor and destructor
  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();

 public:  /// functions
  size_t GetNumSystems() const { return systems_.size(); }
  std::shared_ptr<MicrophysicalSchemeBase> GetSystem(int i) const {
    return systems_[i];
  }
  // void AddFrictionalHeating(std::vector<AirParcel> &air_column) const;

  //! \brief Evolve all microphysical systems
  //!
  //! \param [in,out] ac air column to be evolved
  //! \param [in] time current simulation time
  //! \param [in] dt time step
  void EvolveSystems(AirColumn &ac, Real time, Real dt);

  void Evolve(Real time, Real dt);

  template <typename T>
  void SetConserved(T u, T s);

  template <typename T>
  void GetConserved(T u, T s);

 public:  /// inbound functions
  void SetVsedFromConserved(Hydro const *phydro);

 protected:
  //! sedimentation velocity at cell center [m/s]
  AthenaArray<Real> vsed_[3];

  //! pointers of microphysical systems
  std::vector<std::shared_ptr<MicrophysicalSchemeBase>> systems_;

 private:
  //! meshblock pointer
  MeshBlock const *pmy_block_;

  //! sedimentation options
  SedimentationOptions sed_opts_;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
