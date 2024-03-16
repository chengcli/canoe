#ifndef SRC_HARP_SPECTRAL_GRID_HPP_
#define SRC_HARP_SPECTRAL_GRID_HPP_

// C/C++
#include <memory>
#include <string>
#include <utility>
#include <vector>

// athena
#include <athena/athena.hpp>

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

//! \brief Parse the unit of the spectral grid
//! \return [wavenumber, wavelength, frequency]
std::string parse_unit_with_default(YAML::Node const& my);

//! \brief Read the spectral range from a YAML node
//! \return [wave min, wave max]
std::pair<Real, Real> parse_range(YAML::Node const& my);

//! Base class for a collection of spectral grids
class SpectralGridBase {
 public:
  explicit SpectralGridBase(YAML::Node const& my);

  //! spectral grids
  std::vector<SpectralBin> spec;

 protected:
  //! \brief defines the unit of the spectral grid
  //!
  //! Choose from [wavenumber, wavelength, frequency]
  std::string unit_;
};

//! Policy class for creating a regular spacing spectral grid
class RegularSpacingSpectralGrid : public SpectralGridBase {
 public:
  explicit RegularSpacingSpectralGrid(YAML::Node const& my);
};

//! Policy class for creating a custom spacing spectral grid
class CustomSpacingSpectralGrid : public SpectralGridBase {
 public:
  explicit CustomSpacingSpectralGrid(YAML::Node const& my);
};

//! Policy class for creating a correlated k-table spectral grid
class CKTableSpectralGrid : public SpectralGridBase {
 public:
  explicit CKTableSpectralGrid(YAML::Node const& my);
};

using SpectralGridPtr = std::shared_ptr<SpectralGridBase>;

class SpectralGridFactory {
 public:
  static SpectralGridPtr CreateFrom(YAML::Node const& my);
};

#endif  // SRC_HARP_SPECTRAL_GRID_HPP_
