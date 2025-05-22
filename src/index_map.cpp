// C/C++
#include <mutex>
#include <string>

// canoe
#include <configure.h>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/vectorize.hpp>

// canoe
#include "index_map.hpp"

static std::mutex imap_mutex;

IndexMap::~IndexMap() {
  Application::Logger app("canoe");
  app->Log("Destroy IndexMap");
}

IndexMap const* IndexMap::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(imap_mutex);

  if (myindex_map_ == nullptr) {
    myindex_map_ = new IndexMap();
  }

  return myindex_map_;
}

IndexMap const* IndexMap::InitFromNames(
    std::vector<std::string> const& vapors,
    std::vector<std::string> const& clouds,
    std::vector<std::string> const& tracers) {
  if (myindex_map_ != nullptr) {
    throw RuntimeError("IndexMap", "IndexMap has been initialized");
  }

  myindex_map_ = new IndexMap();

  Application::Logger app("canoe");
  app->Log("Initialize IndexMap");

  // vapor id
  if (vapors.size() > NVAPOR)
    throw ValueError("IndexMap", "Number of vapors", NVAPOR, vapors.size());
  for (size_t i = 0; i < vapors.size(); ++i)
    myindex_map_->vapor_index_map_[vapors[i]] = 1 + i;

  // cloud id
  if (clouds.size() > NCLOUD)
    throw ValueError("IndexMap", "Number of clouds", NCLOUD, clouds.size());
  for (size_t i = 0; i < clouds.size(); ++i)
    myindex_map_->cloud_index_map_[clouds[i]] = i;

  // tracer id
  if (tracers.size() > NTRACER)
    throw ValueError("IndexMap", "Number of tracers", NTRACER, tracers.size());
  for (size_t i = 0; i < tracers.size(); ++i)
    myindex_map_->tracer_index_map_[tracers[i]] = i;

  //! \todo add particle
  return myindex_map_;
}

IndexMap const* IndexMap::InitFromAthenaInput(ParameterInput* pin) {
  if (myindex_map_ != nullptr) {
    throw RuntimeError("IndexMap", "IndexMap has been initialized");
  }

  myindex_map_ = new IndexMap();

  Application::Logger app("canoe");
  app->Log("Initialize IndexMap");

  // vapor id
  std::string str = pin->GetOrAddString("species", "vapor", "");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NVAPOR) {
    throw ValueError("IndexMap", "Number of vapors", NVAPOR, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    myindex_map_->vapor_index_map_[names[i]] = 1 + i;
  }

  // cloud id
  str = pin->GetOrAddString("species", "cloud", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NCLOUD) {
    throw ValueError("IndexMap", "Number of clouds", NCLOUD, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    myindex_map_->cloud_index_map_[names[i]] = i;
  }

  // chemistry id
  str = pin->GetOrAddString("species", "chemistry", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NCHEMISTRY) {
    throw ValueError("IndexMap", "Number of chemistry", NCHEMISTRY,
                     names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    myindex_map_->chemistry_index_map_[names[i]] = i;
  }

  // tracer id
  str = pin->GetOrAddString("species", "tracer", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  if (names.size() > NTRACER) {
    throw ValueError("IndexMap", "Number of tracers", NTRACER, names.size());
  }

  for (size_t i = 0; i < names.size(); ++i) {
    myindex_map_->tracer_index_map_[names[i]] = i;
  }

  // particle id
  str = pin->GetOrAddString("species", "particle", "");
  names = Vectorize<std::string>(str.c_str(), " ,");

  for (size_t i = 0; i < names.size(); ++i) {
    myindex_map_->particle_index_map_[names[i]] = i;
  }

  return myindex_map_;
}

size_t IndexMap::GetSpeciesId(std::string category_name) const {
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
  } else if (category == "chemistry") {
    return NHYDRO + NCLOUD + GetChemistryId(name);
  } else if (category == "tracer") {
    return NHYDRO + NCLOUD + NCHEMISTRY + GetTracerId(name);
  } else {
    throw NotFoundError("GetSpeciesId", "Category " + category);
  }
}

void IndexMap::Destroy() {
  std::unique_lock<std::mutex> lock(imap_mutex);

  if (IndexMap::myindex_map_ != nullptr) {
    delete IndexMap::myindex_map_;
    IndexMap::myindex_map_ = nullptr;
  }
}

std::string IndexMap::GetVaporName(size_t i) const {
  for (auto const& [name, id] : vapor_index_map_) {
    if (id == i) return name;
  }
  throw NotFoundError("GetVaporName", "Vapor id " + std::to_string(i));
}

std::string IndexMap::GetCloudName(size_t i) const {
  for (auto const& [name, id] : cloud_index_map_) {
    if (id == i) return name;
  }
  throw NotFoundError("GetCloudName", "Cloud id " + std::to_string(i));
}

std::string IndexMap::GetTracerName(size_t i) const {
  for (auto const& [name, id] : tracer_index_map_) {
    if (id == i) return name;
  }
  throw NotFoundError("GetTracerName", "Tracer id " + std::to_string(i));
}

IndexMap* IndexMap::myindex_map_ = nullptr;
