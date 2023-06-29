// canoe
#include <configure.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/exceptions.hpp>

// utils
#include <utils/vectorize.hpp>

// canoe
#include "index_map.hpp"

MeshBlock::IndexMap::IndexMap(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  // species id
  std::string str = pin->GetOrAddString("species", "vapor", "");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() != NVAPOR) {
    throw InvalidValueError("IndexMap",
                            "Number of vapors != " + std::to_string(NVAPOR));
  }

  for (size_t i = 0; i < names.size(); ++i) {
    vapor_index_map_[names[i]] = 1 + i;
  }

  str = pin->GetOrAddString("species", "tracer", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() != NTRACER) {
    throw InvalidValueError("IndexMap",
                            "Number of tracers != " + std::to_string(NTRACER));
  }

  for (size_t i = 0; i < names.size(); ++i) {
    tracer_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "cloud", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() != NCLOUDS) {
    throw InvalidValueError("IndexMap",
                            "Number of clouds != " + std::to_string(NCLOUDS));
  }

  for (size_t i = 0; i < names.size(); ++i) {
    cloud_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "chemistry", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() != NCHEMISTRY) {
    throw InvalidValueError(
        "IndexMap", "Number of chemistry != " + std::to_string(NCHEMISTRY));
  }

  for (size_t i = 0; i < names.size(); ++i) {
    chemistry_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "static", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() != NSTATIC) {
    throw InvalidValueError(
        "IndexMap", "Number of static variables != " + std::to_string(NSTATIC));
  }

  for (size_t i = 0; i < names.size(); ++i) {
    static_index_map_[names[i]] = i;
  }

  str = pin->GetOrAddString("species", "particle", "");
  names = Vectorize<std::string>(str.c_str(), " ,");
  for (size_t i = 0; i < names.size(); ++i) {
    particle_index_map_[names[i]] = i;
  }
}

size_t MeshBlock::IndexMap::GetSpeciesId(std::string category_name) const {
  std::string delimiter = ".";

  // Find the position of the delimiter
  size_t delimiter_pos = category_name.find(delimiter);

  if (delimiter_pos == std::string::npos) {
    throw NotFoundError("GetSpeciesId",
                        "Delimiter '" + delimiter + "' in " + category_name);
  }

  // Extract the substrings
  std::string category = category_name.substr(0, delimiter_pos);
  std::string name = category_name.substr(delimiter_pos + delimiter.length());

  if (category == "vapor") {
    return GetVaporId(name);
  } else if (category == "tracer") {
    return NHYDRO + GetTracerId(name);
  } else if (category == "cloud") {
    return NHYDRO + NSCALARS + GetCloudId(name);
  } else if (category == "chemistry") {
    return NHYDRO + NSCALARS + NCLOUDS + GetChemistryId(name);
  } else {
    throw NotFoundError("GetSpeciesId", "Category " + category);
  }
}
