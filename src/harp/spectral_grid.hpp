#ifndef SRC_HARP_SPECTRAL_GRID_HPP_
#define SRC_HARP_SPECTRAL_GRID_HPP_

// C/C++
#include <memory>
#include <string>
#include <utility>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

//! radiation direction
struct Direction {
  Real mu, phi;
};

//! \brief Smallest unit of a spectral grid
struct SpectralBin {
  Real rad, wav1, wav2, wght;
};

//! Base class for a collection of spectral grids
class SpectralGridBase {
 public:
  //! spectral grids
  std::vector<SpectralBin> spec;

  //! \brief defines the unit of the spectral grid
  //!
  //! Choose from [wavenumber, wavelength, frequency]
  std::string unit_type;

  //! \brief Read the spectral range from a YAML node
  //! \return [wave min, wave max]
  std::pair<Real, Real> ReadRangeFrom(YAML::Node const& my);
};

//! Policy class for creating a regular spacing spectral grid
class RegularSpacingSpectralGrid : public SpectralGridBase {
 public:
  RegularSpacingSpectralGrid(YAML::Node const& my);
};

//! Policy class for creating a custom spacing spectral grid
class CustomSpacingSpectralGrid : public SpectralGridBase {
 public:
  CustomSpacingSpectralGrid(YAML::Node const& my);
};

//! \brief Policy class for creating a correlated k-table spectral grid
class CorrelatedKTableSpectralGrid : public SpectralGridBase {
 public:
  CorrelatedKTableSpectralGrid(YAML::Node const& my);
};

using SpectralGridPtr = std::shared_ptr<SpectralGridBase>;

class SpectralGridFactory {
 public:
  static SpectralGridPtr CreateFrom(YAML::Node const& my);
};

#endif  // SRC_HARP_SPECTRAL_GRID_HPP_
