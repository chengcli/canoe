// athena
#include <athena/athena.hpp>

// external
#include <yaml-cpp/yaml.h>

#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <index_map.hpp>

// opacity
#include "Giants/freedman.hpp"
#include "Giants/hydrogen_cia.hpp"
#include "Giants/microwave/mwr_absorbers.hpp"
#include "hitran_absorber.hpp"
#include "nitrogen_cia.hpp"
#include "oxygen_cia.hpp"

AbsorberContainer AbsorberFactory::CreateFrom(
    std::vector<std::string> const& names, std::string band_name,
    YAML::Node const& rad) {
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

AbsorberPtr AbsorberFactory::CreateFrom(YAML::Node const& my,
                                        std::string band_name) {
  if (!my["name"]) {
    throw NotFoundError("AbsorberFactory", "'name' field in absorber");
  }

  if (!my["class"]) {
    throw NotFoundError("AbsorberFactory", "'class' field in absorber " +
                                               my["name"].as<std::string>());
  }

  auto ab = createAbsorberPartial(my["name"].as<std::string>(),
                                  my["class"].as<std::string>());

  if (my["data"]) {
    ab->SetOpacityFile(my["data"].as<std::string>());
  } else {
    ab->SetOpacityFile("kcoeff-" + band_name + ".nc");
  }

  if (my["dependent-species"]) {
    auto species = my["dependent-species"].as<std::vector<std::string>>();
    ab->SetSpeciesIndex(species);
  }

  if (my["model"]) {
    ab->SetModel(my["model"].as<std::string>());
  }

  if (my["parameters"]) {
    ab->SetRealsFrom(my["parameters"]);
  }

  ab->CheckFail();

  return ab;
}

namespace gp = GiantPlanets;

AbsorberPtr AbsorberFactory::createAbsorberPartial(std::string name,
                                                   std::string type) {
  AbsorberPtr ab;

  if (type == "XIZ-H2-H2-CIA") {
    ab = std::make_shared<XizH2H2CIA>();
  } else if (type == "XIZ-H2-He-CIA") {
    ab = std::make_shared<XizH2HeCIA>();
  } else if (type == "Hitran") {
    ab = std::make_shared<HitranAbsorber>(name);
  } else if (type == "freedman_simple") {
    ab = std::make_shared<FreedmanSimple>();
  } else if (type == "freedman_mean") {
    ab = std::make_shared<FreedmanSimple>();
  } else if (type == "radio-NH3") {
    ab = std::make_shared<gp::MwrAbsorberNH3>();
  } else if (type == "radio-H2O") {
    ab = std::make_shared<gp::MwrAbsorberH2O>();
  } else if (type == "radio-H2S") {
    ab = std::make_shared<gp::MwrAbsorberH2S>();
  } else if (type == "radio-PH3") {
    ab = std::make_shared<gp::MwrAbsorberPH3>();
  } else if (type == "radio-CIA") {
    ab = std::make_shared<gp::MwrAbsorberCIA>();
  } else if (type == "radio-Electron") {
    ab = std::make_shared<gp::MwrAbsorberElectron>();
  } else {
    throw NotFoundError("createAbsorberPartial", "Class = " + type);
  }

  return ab;
}
