// C/C++
#include <memory>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/exceptions.hpp>

// harp
#include "spectral_grid.hpp"

SpectralGridPtr SpectralGridFactory::CreateFrom(YAML::Node const& my) {
  if (!my["grid-type"]) {
    throw NotFoundError("SpectralGridFactory", "grid-type");
  }

  std::string grid_type = my["grid-type"].as<std::string>();

  SpectralGridPtr pgrid;

  if (grid_type == "regular") {
    pgrid = std::make_shared<RegularSpacingSpectralGrid>(my);
  } else if (grid_type == "custom") {
    pgrid = std::make_shared<CustomSpacingSpectralGrid>(my);
  } else if (grid_type == "cktable") {
    pgrid = std::make_shared<CKTableSpectralGrid>(my);
  } else {
    throw RuntimeError("SpectralGridFactory", "unknown grid type");
  }

  return pgrid;
}
