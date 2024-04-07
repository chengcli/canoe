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

  auto air = AirParcel::gather_from_primitive(pmb, k, ju_, pmb->is);
  air.ToMoleFraction();

  int iter = 0, max_iter = 200;
  while (iter++ < max_iter) {
    // read in vapors
    for (auto n : GetSpeciesIndex()) {
      air.w[n] = par[n];
    }

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate(&air, -dlnp / 2., "dry");
      if (air.w[IPR] < P0) break;
    }

    // extrapolate down to where air is
    pthermo->Extrapolate(&air, log(P0 / air.w[IPR]), "dry");

    // make up for the difference
    Ts += T0 - air.w[IDN];
    if (std::abs(T0 - air.w[IDN]) < 0.01) break;

    // app->Log("Iteration #", iter);
    // app->Log("T", air.w[IDN]);
  }
}
