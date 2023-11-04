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
#include <common.hpp>
#include <virtual_groups.hpp>

// harp
#include "spectral_grid.hpp"

// exchanger
#include <exchanger/linear_exchanger.hpp>

class OutputParameters;
class Absorber;
class SpectralGridBase;
struct Direction;

class RadiationBand : public NamedGroup,
                      public FlagGroup,
                      public ParameterGroup,
                      public ASCIIOutputGroup,
                      public LinearExchanger<RadiationBand> {
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
  RadiationBand(std::string name, YAML::Node const &rad);
  virtual ~RadiationBand();

 public:  // member functions
  //! \brief Allocate memory for radiation band
  void Resize(int nc1, int nc2 = 1, int nc3 = 1);

  //! \brief Create radiative transfer solver from YAML node
  //!
  //! \param[in] my YAML node containing the current band
  std::shared_ptr<RTSolver> CreateRTSolverFrom(YAML::Node const &my);

  //! Get number of spectral grids
  size_t GetNumSpecGrids() const { return pgrid_->spec.size(); }

  //! Get number of absorbers
  size_t GetNumAbsorbers() const { return absorbers_.size(); }

  //! Get number of phase function moments
  size_t GetNumPhaseMoments() const { return nphase_moments; }

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

 public:  // inbound functions
  //! \brief Set spectral properties for an air column
  //!
  //! \param[in] air mole fraction representation of air column
  void SetSpectralProperties(AirColumn &air, Coordinates const *pcoord, int k,
                             int j);

 public:  // ASCIIOutputGroup functions
  void WriteAsciiHeader(OutputParameters const *) const override;
  void WriteAsciiData(OutputParameters const *) const override;

 public:  // ColumnExchanger functions
  void PackData() override;
  bool UnpackData() override;
  void ClearBoundary() override;

 protected:
  //! \brief number of phase function moments
  size_t nphase_moments_;

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
using RadiationBandContainer = std::vector<RadiationBandPtr>;

// Specialization for RadiationBand Exchanger
template <>
struct MessageTraits<RadiationBand> {
  using DataType = Real;

  constexpr static int num_buffers = 2;
  constexpr static std::string name = "RadiationBand";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif  // MPI_PARALLEL
};

class RadiationBandsFactory {
 public:
  static RadiationBandContainer CreateFrom(ParameterInput *pin,
                                           std::string key);
  static RadiationBandContainer CreateFrom(std::string filename);
};

#endif  // SRC_HARP_RADIATION_BAND_HPP_
