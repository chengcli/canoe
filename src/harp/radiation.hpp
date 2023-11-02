#ifndef SRC_HARP_RADIATION_HPP_
#define SRC_HARP_RADIATION_HPP_

// C/C++ headers
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

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

  //! \brief Get total number of incoming rays
  size_t GetNumOutgoingRays() const;

 public:  // inbound functions
  //! \brief Calculate the radiative flux
  void CalRadiativeFlux(MeshBlock const *pmb, Real time, int k, int j, int il,
                        int iu);

  //! \brief Calculate the radiance
  void CalRadiance(MeshBlock const *pmb, Real time, int k, int j, int il,
                   int iu);

 public:  // outbound functions
  //! \brief Add the radiative flux to hydro energy flux
  void AddRadiativeFlux(Hydro *phydro, int k, int j, int il, int iu) const;

 public:  // RestartGroup functions
  size_t RestartDataSizeInBytes() const override;
  size_t DumpRestartData(char *pdst) const override;
  size_t LoadRestartData(char *psrt) override;

 protected:
  //! all radiation bands
  std::vector<RadiationBandPtr> bands_;

  //! incomming rays
  std::vector<Direction> rayInput_;
};

using RadiationPtr = std::shared_ptr<Radiation>;

namespace RadiationFlags {

const uint64_t None = 0LL;
const uint64_t Dynamic = 1LL << 0;
const uint64_t LineByLine = 1LL << 1;
const uint64_t CorrelatedK = 1LL << 2;
const uint64_t Planck = 1LL << 3;
const uint64_t Star = 1LL << 4;
const uint64_t Sphere = 1LL << 5;
const uint64_t FluxOnly = 1LL << 6;
const uint64_t Normalize = 1LL << 7;
const uint64_t WriteBinRadiance = 1LL << 8;

//! \todo Do we put scattering model flag here?
//! e.g. HENYEY_GREENSTEIN, ISOTROPIC, RAYLEIGH, MIE, etc.

}  // namespace RadiationFlags

namespace RadiationHelper {
std::vector<Direction> parse_radiation_directions(std::string str);
uint64_t parse_radiation_flags(std::string str);
void get_phase_momentum(Real *pmom, int iphas, Real gg, int npmom);
};  // namespace RadiationHelper

#endif  //  SRC_HARP_RADIATION_HPP_
