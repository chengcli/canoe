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
  // vapor id
  std::string str = pin->GetOrAddString("species", "vapor", "");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NVAPOR) {
    throw RuntimeError("IndexMap", "Number of vapors", NVAPOR, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    vapor_index_map_[names[i]] = 1 + i;
  }

  // cloud id
  str = pin->GetOrAddString("species", "cloud", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NCLOUD) {
    throw RuntimeError("IndexMap", "Number of clouds", NCLOUD, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    cloud_index_map_[names[i]] = i;
  }

  // tracer id
  str = pin->GetOrAddString("species", "tracer", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NTRACER) {
    throw RuntimeError("IndexMap", "Number of tracers", NTRACER, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    tracer_index_map_[names[i]] = i;
  }

  // chemistry id
  str = pin->GetOrAddString("species", "chemistry", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NCHEMISTRY) {
    throw RuntimeError("IndexMap", "Number of chemistry", NCHEMISTRY,
                       names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    chemistry_index_map_[names[i]] = i;
  }

  // static variable id
  str = pin->GetOrAddString("species", "static", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NSTATIC) {
    throw RuntimeError("IndexMap", "Number of static variables", NSTATIC,
                       names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    static_index_map_[names[i]] = i;
  }

  // particle id
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
  } else if (category == "cloud") {
    return NHYDRO + GetCloudId(name);
  } else if (category == "tracer") {
    return NHYDRO + NCLOUD + GetTracerId(name);
  } else if (category == "chemistry") {
    return NHYDRO + NCLOUD + NTRACER + GetChemistryId(name);
  } else {
    throw NotFoundError("GetSpeciesId", "Category " + category);
  }
}
