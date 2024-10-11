#pragma once

// C/C++
#include <map>
#include <memory>
#include <vector>

// athena
#include <athena/athena.hpp>

// snap
#include <snap/stride_iterator.hpp>

class MeshBlock;
class ParameterInput;

namespace Cantera {
class Kinetics;
}  // namespace Cantera

//! \brief root-level management class for microphysics
class Microphysics {
 public:  /// public access members
  friend std::vector<std::shared_ptr<Cantera::Kinetics>> get_kinetics_objects(
      Microphysics const *pmicro);

  //! \todo(CLI) track cloud temperature and momentum
  //! tem, v1, v2, v3

  //! microphysics input key in the input file [microphysics_config]
  static const std::string input_key;

  //! sedimentation velocity at cell interface [m/s]
  AthenaArray<Real> vsedf[3];

 public:  /// constructor and destructor
  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();

 public:  /// functions
  // void AddFrictionalHeating(std::vector<AirParcel> &air_column) const;

  //! \brief Evolve all microphysical systems
  //!
  //! \param [in] time current simulation time
  //! \param [in] dt time step
  void Evolve(Real time, Real dt);

  template <typename T>
  void SetConserved(T u, T s, T m);

  template <typename T>
  void GetConserved(T u, T s, T m);

 public:
  template <typename T>
  void ModifyScalarFlux(T sflux[3]);

 protected:
  //! sedimentation velocity at cell center [m/s]
  AthenaArray<Real> vsed_[3];

  //! pointers of microphysical systems
  std::vector<std::shared_ptr<Cantera::Kinetics>> systems_;

 private:
  //! meshblock pointer
  MeshBlock const *pmy_block_;
};

using MicrophysicsPtr = std::shared_ptr<Microphysics>;
