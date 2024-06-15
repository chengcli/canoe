// C/C++
#include <memory>
#include <string>

// cantera
#include <cantera/base/yaml.h>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// canoe
#include <air_parcel.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

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
  auto pthermo = Thermodynamics::GetInstance();

  auto air = AirParcelHelper::gather_from_primitive(pmb, k, ju_, pmb->is);
  air.ToMoleFraction();

  Real dlnp = 1.;

  // read in vapors
  for (auto n : GetMySpeciesIndices()) {
    air.w[n] = par[n];
  }

  for (int i = pmb->is; i <= pmb->ie; ++i) {
    pthermo->Extrapolate(&air, -dlnp / 2., "dry");
  }
}
