// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/athena.hpp>

// canoe
#include <air_parcel.hpp>
#include <index_map.hpp>

// opacity
#include "absorber.hpp"

Absorber::Absorber(std::string name, std::vector<std::string> const& names)
    : NamedGroup(name) {
  Application::Logger app("harp");
  app->Log("Create Absorber " + name_);

  auto pindex = IndexMap::GetInstance();

  for (auto& s : names) {
    imols_.push_back(pindex->GetSpeciesId(s));
  }

  app->Log("Dependent species ids", imols_);
}

Absorber::~Absorber() {
  Application::Logger app("harp");
  app->Log("Destroy Absorber " + name_);
}

AbsorberContainer AbsorberFactory::CreateFrom(
    std::vector<std::string> const& names, std::string band_name,
    YAML::Node& rad) {
  AbsorberContainer absorbers;

  for (auto& name : names) {
    bool found = false;
    for (auto& my : rad["opacity-sources"]) {
      if (name == my["name"].as<std::string>()) {
        absorbers.push_back(AbsorberFactory::CreateFrom(my, band_name));
        found = true;
        break;
      }
    }

    if (!found) {
      throw NotFoundError("AbsorberFactory", "Opacity " + name);
    }
  }

  return absorbers;
}

AbsorberPtr AbsorberFactory::CreateFrom(YAML::Node& my, std::string band_name) {
  std::vector<std::string> species;
  std::string name = my["name"].as<std::string>();

  std::string data_file;
  std::string full_path;

  if (my["dependent-species"]) {
    species = my["dependent-species"].as<std::vector<std::string>>();
  }

  if (my["data"]) {
    data_file = my["data"].as<std::string>();
  } else {
    data_file = "kcoeff-" + band_name + ".nc";
  }

  auto app = Application::GetInstance();

  try {
    full_path = app->FindInputFile(data_file);
  } catch (NotFoundError const& e) {
    auto log = app->GetMonitor("harp");
    std::stringstream ss;
    ss << e.what() << std::endl;
    ss << name << " will not be loaded." << std::endl;
    log->Warn(ss.str());
  }

  auto ab = CreateAbsorber(name, species, full_path);

  if (my["model"]) {
    ab->SetModel(my["model"].as<std::string>());
  }

  if (my["parameters"]) {
    ab->SetRealsFrom(my["parameter"]);
  }

  return ab;
}

void __attribute__((weak))
AbsorberFactory::CreateAbsorber(std::string name,
                                std::vector<std::string> const& species,
                                std::string full_path) {}
