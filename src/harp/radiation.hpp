#ifndef SRC_HARP_RADIATION_HPP_
#define SRC_HARP_RADIATION_HPP_

// C/C++ headers
#include <memory>
#include <string>
#include <vector>

// astro
#include <astro/celestrial_body.hpp>

// athena
#include <athena/athena.hpp>

// harp
#include "radiation_band.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;

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
}  // namespace RadiationFlags

class Radiation {
 public:
  // access data
  AthenaArray<Real> radiance;

  // radiation bands
  std::vector<RadiationBandPtr> bands;

  AthenaArray<Real> flxup, flxdn;

  // functions
  Radiation(MeshBlock *pmb, ParameterInput *pin);

  ~Radiation();

  size_t GetNumBands() { return bands.size(); }

  RadiationBand *GetBand(int i) { return bands[i].get(); }

  RadiationBand *GetBand(std::string const &name);

  void PopulateRadiationBands(ParameterInput *pin);

  void CalRadiativeFlux(Real time, int k, int j, int il, int iu);

  void CalRadiance(Real time, int k, int j, int il, int iu);

  void addRadiativeFlux(Hydro *phydro, int k, int j, int il, int iu) const;

  size_t GetNumOutgoingRays() const;

  int test(uint64_t flag) const { return rflags_ & flag; }

  void set(uint64_t flag) { rflags_ |= flag; }

  // restart functions
  size_t getRestartDataSizeInBytes() const;

  size_t dumpRestartData(char *pdst) const;

  size_t loadRestartData(char *psrc);

 protected:
  // data
  uint64_t rflags_;
  Real cooldown_, current_;

  // incomming rays
  std::vector<Direction> rayInput_;

  Real stellarDistance_au_;

  // connections
  MeshBlock *pmy_block_;
  CelestrialBodyPtr planet_;
};

using RadiationPtr = std::shared_ptr<Radiation>;

#endif  //  SRC_HARP_RADIATION_HPP_
