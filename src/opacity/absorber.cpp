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
    std::vector<std::string> const& names, YAML::Node& rad) {
  AbsorberContainer absorbers;

  for (auto& name : names) {
    bool found = false;
    for (auto& my : rad["opacity-sources"]) {
      if (name == my["name"].as<std::string>()) {
        absorbers.push_back(AbsorberFactory::CreateFrom(my));
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

AbsorberPtr AbsorberFactory::CreateFrom(YAML::Node& my) {
  AbsorberPtr ab;

  if (!my["parameters"]) {
    my["parameters"] = {};
  }

  if (!my["model"]) {
    my["model"] = "unspecified";
  }

  ab->SetModel(my["model"].as<std::string>());

  return ab;
}
