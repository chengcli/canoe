#ifndef SRC_HARP_RADIATION_HPP_
#define SRC_HARP_RADIATION_HPP_

// C/C++ headers
#include <memory>
#include <string>
#include <string_view>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <common.hpp>  // Direction

// harp
#include "radiation_band.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;
class RadiationBand;

class Radiation : public RestartGroup,
                  // public MeshOutputGroup,
                  public CounterGroup,
                  public FlagGroup {
 public:  // public access data
  //! radiation input key in the input file [radiation_config]
  static const std::string input_key;

  //! radiance of all bands
  AthenaArray<Real> radiance;

  //! upward flux of all bands
  AthenaArray<Real> flxup;

  //! downward flux of all bands
  AthenaArray<Real> flxdn;

 public:  // constructor and destructor
  // Radiation() {}
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

 public:  // member functions
  //! Get number of bands
  size_t GetNumBands() const { return bands_.size(); }

  //! Get band by index
  std::shared_ptr<RadiationBand> GetBand(int i) const { return bands_[i]; }

  //! Get band by name
  std::shared_ptr<RadiationBand> GetBandByName(std::string const &name) const;

  //! Get incoming ray by index
  Direction const &GetRayInput(int i) const { return rayInput_[i]; }

  //! \brief Get total number of incoming rays
  size_t GetNumOutgoingRays() const;

  //! Get the default surface temperature (K)
  Real GetSurfaceTemperature() const { return surface_temperature_; }

 public:  // inbound functions
  //! \brief Calculate the radiative flux
  void CalFlux(MeshBlock const *pmb, int k, int j, int il, int iu);

  //! \brief Calculate the radiance
  void CalRadiance(MeshBlock const *pmb, int k, int j);

 public:  // outbound functions
  //! \brief Add the radiative flux to hydro energy flux
  void AddRadiativeFlux(Hydro *phydro, int k, int j, int il, int iu) const;

 public:  // RestartGroup functions
  size_t RestartDataSizeInBytes(Mesh const *pm) const override;
  void DumpRestartData(char *pdst) const override;
  size_t LoadRestartData(char *psrt) override;

 protected:
  //! all radiation bands
  std::vector<RadiationBandPtr> bands_;

  //! incomming rays
  std::vector<Direction> rayInput_;

  //! surface temperature
  Real surface_temperature_;
};

using RadiationPtr = std::shared_ptr<Radiation>;

namespace RadiationFlags {

const uint64_t None = 0LL;
const uint64_t TimeDependent = 1LL << 0;
const uint64_t BroadBand = 1LL << 1;
const uint64_t ThermalEmission = 1LL << 2;
const uint64_t StellarBeam = 1LL << 3;
const uint64_t Normalize = 1LL << 7;
const uint64_t WriteBinRadiance = 1LL << 8;

}  // namespace RadiationFlags

namespace RadiationHelper {
//! \brief Get the number of grids in the outgoing ray directions
std::pair<std::vector<Real>, std::vector<Real>> get_direction_grids(
    std::vector<Direction> const &dirs);

//! \brief Parse radiation direction string
//!
//! First number is the polar angle (degrees), second number is the azimuthal
//! angle (degrees) \param[in] str radiation direction string, e.g. (45, 30)
//! \return radiation direction
Direction parse_radiation_direction(std::string_view str);

//! \brief Parse radiation directions string, sperated by comma
//!
//! Example input string: "(45, 30), (45, 60)"
//! \param[in] str radiation directions string
std::vector<Direction> parse_radiation_directions(std::string str);
uint64_t parse_radiation_flags(std::string str);
void get_phase_momentum(Real *pmom, int iphas, Real gg, int npmom);
};  // namespace RadiationHelper

#endif  //  SRC_HARP_RADIATION_HPP_
