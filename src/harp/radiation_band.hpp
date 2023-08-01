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
  AthenaArray<Real> bflxup, bflxdn;

  // btoa is a reference to radiance in Radiation
  AthenaArray<Real> btoa;

  // functions
  RadiationBand() {}

  RadiationBand(MeshBlock *pmb, ParameterInput *pin, YAML::Node &node,
                std::string name);

  ~RadiationBand();

  size_t GetNumBins() const { return spec_.size(); }

  Real GetWavenumberMin() const { return wmin_; }

  Real GetWavenumberMax() const { return wmax_; }

  Real GetWavenumberRes() const;

  bool HasParameter(std::string const &name) const {
    return params_.count(name);
  }

  Real GetParameter(std::string const &name) const { return params_.at(name); }

  size_t GetNumAbsorbers() const { return absorbers_.size(); }

  AbsorberPtr GetAbsorber(int i) { return absorbers_[i]; }

  AbsorberPtr GetAbsorberByName(std::string const &name);

  size_t GetNumOutgoingRays() { return rayOutput_.size(); }

  std::string GetName() { return name_; }

  std::string GetCategory() { return category_; }

  Real GetCosinePolarAngle(int n) const { return rayOutput_[n].mu; }

  Real GetAzimuthalAngle(int n) const { return rayOutput_[n].phi; }

  void AddAbsorber(ParameterInput *pin, YAML::Node &node);

  void SetSpectralProperties(int k, int j, int il, int iu);

  void WriteBinRadiance(OutputParameters const *) const;

  // implementation of RT Solver
  class RTSolver;
  class RTSolverLambert;
  class RTSolverDisort;
  std::shared_ptr<RTSolver> psolver;

 protected:
  void setWavenumberRange(YAML::Node &my);
  void setWavenumberGrid(YAML::Node &my);
  void setFrequencyRange(YAML::Node &my);
  void setFrequencyGrid(YAML::Node &my);
  void setWavelengthRange(YAML::Node &my);
  void setWavelengthGrid(YAML::Node &my);

  void addAbsorberGiants(ParameterInput *pin, YAML::Node &node);
  void addAbsorberGiantsRadio(YAML::Node &node);
  void addAbsorberGiantsInfrared(YAML::Node &node);
  void addAbsorberGiantsVisible(YAML::Node &node);
  void addAbsorberGiantsUltraviolet(YAML::Node &node);

  void addAbsorberEarth(ParameterInput *pin, YAML::Node &node);
  void addAbsorberEarthInfrared(YAML::Node &node);
  void addAbsorberEarthVisible(YAML::Node &node);
  void addAbsorberEarthUltraviolet(YAML::Node &node);

  void addAbsorberVenus(ParameterInput *pin, YAML::Node &node);

  void addAbsorberMars(ParameterInput *pin, YAML::Node &node);
  void addAbsorberMarsInfrared(YAML::Node &node);

  int test(uint64_t flag) const { return bflags_ & flag; }
  void set(uint64_t flag) { bflags_ |= flag; }

  // absorbers
  std::vector<AbsorberPtr> absorbers_;

  // data
  std::string name_;
  std::string myfile_;
  std::string category_;
  uint64_t bflags_;

  // spectra
  Real wmin_, wmax_;
  std::vector<Spectrum> spec_;

  // outgoing rays
  std::vector<Direction> rayOutput_;

  // bin data
  // 2d
  AthenaArray<Real> tau_, ssa_, toa_;
  AthenaArray<Real> flxup_, flxdn_;
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
