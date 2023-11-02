#ifndef SRC_HARP_SPECTRAL_GRID_HPP_
#define SRC_HARP_SPECTRAL_GRID_HPP_

// C/C++
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
  std::pair<Real, Real> ReadRangeFrom(YAML::Node &node);

  //! number of the spectral bin
  size_t GetSize() const { return spec.size(); }
};

//! Policy class for creating a regular spacing spectral grid
class RegularSpacingSpectralGrid : public SpectralGridBase {
 public:
  RegularSpacingSpectralGrid(YAML::Node &my);
};

//! Policy class for creating a custom spacing spectral grid
class CustomSpacingSpectralGrid : public SpectralGridBase {
 public:
  CustomSpacingSpectralGrid(YAML::Node &my);
};

//! \brief Policy class for creating a correlated k-table spectral grid
class CorrelatedKTableSpectralGrid : public SpectralGridBase {
 public:
  CorrelatedKTableSpectralGrid(YAML::Node &my);
};

using SpectralGridPtr = std::shared_ptr<SpectralGridBase>;

class SpectralGridFactory {
 public:
  static SpectralGridPtr = CreateFrom(YAML::Node & my);
};

#endif  // SRC_HARP_SPECTRAL_GRID_HPP_
