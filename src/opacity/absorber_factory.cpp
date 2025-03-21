// athena
#include <athena/athena.hpp>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <index_map.hpp>

// opacity
#include "Giants/grey_gas.hpp"
#include "Giants/hydrogen_cia.hpp"
#include "Giants/microwave/mwr_absorbers.hpp"
#include "absorber_ck.hpp"
#include "hitran_absorber.hpp"
#include "nitrogen_cia.hpp"
#include "oxygen_cia.hpp"
#include "simple_cloud.hpp"

AbsorberContainer AbsorberFactory::CreateFrom(
    std::vector<std::string> const& names, std::string band_name,
    YAML::Node const& rad) {
  AbsorberContainer absorbers;

  for (auto& name : names) {
    AbsorberContainer atmp;

    for (auto& my : rad["opacity-sources"]) {
      if (name == my["name"].as<std::string>()) {
        atmp.push_back(AbsorberFactory::CreateFrom(my, band_name));
      }
    }

    if (atmp.empty()) {
      throw NotFoundError("AbsorberFactory", "Opacity " + name);
    } else if (atmp.size() >
               1) {  // use band_name to uniquely identify the absorber
      bool found = false;
      for (auto& ab : atmp) {
        if (ab->GetOpacityFile().find(band_name) != std::string::npos) {
          absorbers.push_back(ab);
          found = true;
          break;
        }
      }

      if (!found) {
        throw NotFoundError("AbsorberFactory",
                            "There are multiple absorbers with the same name " +
                                name +
                                ". Rename the `data` field in the absorber to "
                                "identify the band.");
      }
    } else {
      absorbers.push_back(atmp[0]);
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
    ab->SetOpacityFile("");
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
  } else if (type == "HitranCK") {
    ab = std::make_shared<HitranAbsorberCK>(name);
  } else if (type == "SimpleGrey") {
    ab = std::make_shared<SimpleGrey>(name);
  } else if (type == "FreedmanMean") {
    ab = std::make_shared<FreedmanMean>(name);
  } else if (type == "JupGasv") {
    ab = std::make_shared<JupGasv>(name);
  } else if (type == "JupGasir") {
    ab = std::make_shared<JupGasir>(name);
  } else if (type == "HeliosCK") {
    ab = std::make_shared<HeliosCK>(name);
  } else if (type == "SimpleCloud") {
    ab = std::make_shared<SimpleCloud>(name);
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
