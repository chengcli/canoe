// C/C++
#include <memory>
#include <string>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"

CompositionInversion::CompositionInversion(YAML::Node const &node)
    : Inversion("composition") {
  Application::Logger app("inversion");
  app->Log("Initializing CompositionInversion");
  char buf[80];

  // species id
  auto species = node["variables"].as<std::vector<std::string>>();
  SetSpeciesIndex(species);
}

void CompositionInversion::UpdateModel(MeshBlock *pmb,
                                       std::vector<Real> const &par,
                                       int k) const {
  Application::Logger app("inversion");
  app->Log("Update Model");

  // int is = pblock_->is, ie = pblock_->ie;
}
