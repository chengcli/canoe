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
#include <athena/athena.hpp>  // Real

// canoe
#include <air_parcel.hpp>
#include <virtual_groups.hpp>

class OutputParameters;
class Absorber;
struct Direction;

class RadiationBand : public NamedGroup,
                      public FlagGroup,
                      public ParameterGroup,
                      public ASCIIOutputGroup {
 public:  // public access data
  // implementation of RT Solver
  class RTSolver;
  class RTSolverLambert;
  class RTSolverDisort;

  //! radiative transfer solver
  std::shared_ptr<RTSolver> psolver;

  //! band optical depth
  AthenaArray<Real> btau;

  //! band single scattering albedo
  AthenaArray<Real> bssa;

  //! band phase function moments
  AthenaArray<Real> bpmom;

  //! band upward flux
  AthenaArray<Real> bflxup;

  //! band downward flux
  AthenaArray<Real> bflxdn;

  //! \brief band top-of-the-atmosphere radiance
  //!
  //! This is a shallow reference to radiance in Radiation
  AthenaArray<Real> btoa;

 public:  // constructor and destructor
  RadiationBand() {}
  RadiationBand(YAML::Node &inp, std::string name, int nc1 = 1, int nc2 = 1,
                int nc3 = 1);
  ~RadiationBand();

 public:  // functions
  //! \brief Create radiative transfer solver from YAML node
  //!
  //! \param[in] my YAML node containing the current band
  static std::shared_ptr<RTSolver> CreateRTSolverFrom(YAML::Node &my);

  //! Get number of spectral grids
  size_t GetNumSpecGrids() const { return spec_grid_->GetSize(); }

  //! Get number of absorbers
  size_t GetNumAbsorbers() const { return absorbers_.size(); }

  //! Get an individual absorber
  std::shared_ptr<Absorber> GetAbsorber(int i) { return absorbers_[i]; }

  //! Get an individual absorber by name
  std::shared_ptr<Absorber> GetAbsorberByName(std::string const &name);

  //! Get number of outgoing rays
  size_t GetNumOutgoingRays() { return rayOutput_.size(); }

  //! Get cosine polar angle of the outgoing ray
  Real GetCosinePolarAngle(int n) const { return rayOutput_[n].mu; }

  //! Get azimuthal angle of the outgoing ray
  Real GetAzimuthalAngle(int n) const { return rayOutput_[n].phi; }

 public:  // incoming functions
  //! \brief Set spectral properties for an air column
  //!
  //! \param[in] air mole fraction representation of air column
  void SetSpectralProperties(AirColumn const &air, Coordinates const *pcoord,
                             int k, int j);

  void WriteBinRadiance(OutputParameters const *) const;

 protected:
  //! spectral grid
  std::shared_ptr<SpectralGridBase> pgrid_;

  //! \brief identify wave- number (wavelength, frequency) range
  std::pair<Real, Real> wrange_;

  //! all absorbers
  std::vector<std::shared_ptr<Absorber>> absorbers_;

  //! outgoing rays
  std::vector<Direction> rayOutput_;

  //! \brief spectral bin optical depth
  //!
  //! This is a two-dimensional array. The first dimension is the size of the
  //! spectral bin and the second dimension is the size of the vertical
  //! dimension.
  AthenaArray<Real> tau_;

  //! \brief spectral bin single scattering albedo
  //!
  //! This is a two-dimensional array. The first dimension is the size of the
  //! spectral bin and the second dimension is the size of the vertical
  //! dimension.
  AthenaArray<Real> ssa_;

  //! \brief spectral bin top-of-the-atmosphere radiance
  //!
  //! This is a two-dimensional array. The first dimension is the size of the
  //! spectral bin and the second dimension is the size of the vertical
  //! dimension.
  AthenaArray<Real> toa_;

  //! \brief spectral bin upward flux
  //!
  //! This is a two-dimensional array. The first dimension is the size of the
  //! spectral bin and the second dimension is the size of the vertical
  //! dimension.
  AthenaArray<Real> flxup_;

  //! \brief spectral bin downward flux
  //!
  //! This is a two-dimensional array. The first dimension is the size of the
  //! spectral bin and the second dimension is the size of the vertical
  //! dimension.
  AthenaArray<Real> flxdn_;

  //! \brief spectral bin phase function moments
  //!
  //! This is a three-dimensional array. The first dimension is the size of the
  //! spectral bin. The second dimension is the size of the vertical dimension.
  //! The third dimension is the size of the phase function moments (1 + npmom).
  AthenaArray<Real> pmom_;

  //! temperature at cell center
  AthenaArray<Real> tem_;

  //! temperature at cell boundary (face)
  AthenaArray<Real> temf_;
};

using RadiationBandPtr = std::shared_ptr<RadiationBand>;

class RadiationBandsFactory {
 public:
  static std::vector<RadiationBandPtr> CreateFrom(ParameterInput *pin,
                                                  std::string key);

  static std::vector<RadiationBandPtr> CreateFrom(std::string filename);
};

#endif  // SRC_HARP_RADIATION_BAND_HPP_
