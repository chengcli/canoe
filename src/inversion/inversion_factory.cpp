// C/C++
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// utils
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"

std::vector<InversionPtr> InversionFactory::CreateFrom(YAML::Node const& node) {
  Application::Logger app("inversion");
  app->Log("Create inversion queue");

  if (!node["tasks"]) {
    throw NotFoundError("InversionFactory::CreateFrom", "tasks");
  }

  if (!node["tasks"].IsSequence()) {
    throw RuntimeError("InversionFactory::CreateFrom",
                       "tasks must be a sequence");
  }

  std::vector<std::string> task_names;

  for (auto const& item : node["tasks"]) {
    task_names.push_back(item.as<std::string>());
  }

  std::vector<InversionPtr> all_fits;
  InversionPtr pfit;

  int jl = NGHOST, ju = NGHOST;
  for (auto p : task_names) {
    if (p == "profile") {
      pfit = std::make_shared<ProfileInversion>(node["profile"]);
    } else if (p == "composition") {
      pfit = std::make_shared<CompositionInversion>(node["composition"]);
    } else {
      throw NotFoundError("CreateAllInversions", p);
    }

    ju += pfit->GetSteps();
    pfit->SetStepRange(jl, ju);
    jl = ju;
    all_fits.push_back(std::move(pfit));
  }

  app->Log("Number of inversions = " + std::to_string(all_fits.size()));

  return all_fits;
}
