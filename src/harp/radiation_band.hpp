#ifndef SRC_HARP_RADIATION_BAND_HPP_
#define SRC_HARP_RADIATION_BAND_HPP_

// C/C++ headers
#include <map>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>

// harp
#include "absorber.hpp"

struct Spectrum {
  Real rad, wav1, wav2, wght;
};

struct Direction {
  Real mu, phi;
};

class OutputParameters;

class RadiationBand {
 public:
  // access data
  AthenaArray<Real> btau, bssa, bpmom;

  // functions
  RadiationBand(MeshBlock *pmb, ParameterInput *pin, YAML::Node &node,
                std::string name);

  ~RadiationBand();

  size_t GetNumAbsorbers() { return absorbers_.size(); }

  Absorber *GetAbsorber(int i) { return absorbers_[i].get(); }

  Absorber *GetAbsorber(std::string const &name);

  size_t getNumOutgoingRays() { return rayOutput_.size(); }

  std::string GetName() { return name_; }

  Real getCosinePolarAngle(int n) const { return rayOutput_[n].mu; }

  Real getAzimuthalAngle(int n) const { return rayOutput_[n].phi; }

  void AddAbsorber(ParameterInput *pin, std::string bname, YAML::Node &node);

  void SetSpectralProperties(int k, int j, int il, int iu);

  void calculateBandFlux(AthenaArray<Real> *flxup, AthenaArray<Real> *flxdn,
                         Direction const &rayInput, Real dist_au, int k, int j,
                         int il, int iu);

  void calculateBandRadiance(AthenaArray<Real> *radiance,
                             Direction const &rayInput, Real dist_au, int k,
                             int j, int il, int iu);

  int test(uint64_t flag) const { return bflags_ & flag; }

  void set(uint64_t flag) { bflags_ |= flag; }

  void writeBinRadiance(OutputParameters const *) const;

  // user implementation of RT Solver
  class Impl;
  std::shared_ptr<Impl> pimpl;

 protected:
  void addAbsorberGiants(ParameterInput *pin, std::string bname,
                         YAML::Node &node);
  void addAbsorberEarth(ParameterInput *pin, std::string bname,
                        YAML::Node &node);
  void addAbsorberVenus(ParameterInput *pin, std::string bname,
                        YAML::Node &node);
  void addAbsorberMars(ParameterInput *pin, std::string bname,
                       YAML::Node &node);

  // absorbers
  std::vector<AbsorberPtr> absorbers_;

  // data
  std::string name_;
  uint64_t bflags_;

  // spectra
  Real wmin_, wmax_;
  std::vector<Spectrum> spec_;

  // outgoing rays
  std::vector<Direction> rayOutput_;

  // bin data
  // 2d
  AthenaArray<Real> tau_, ssa_, toa_;
  // 3d
  AthenaArray<Real> pmom_;
  // 1d
  AthenaArray<Real> tem_, temf_;

  Real alpha_;  // T ~ Ts*(\tau/\tau_s)^\alpha at lower boundary
  ParameterMap params_;

  // connection
  MeshBlock *pmy_block_;
};

using RadiationBandPtr = std::shared_ptr<RadiationBand>;

#endif  // SRC_HARP_RADIATION_BAND_HPP_
