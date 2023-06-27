// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/cell_variables.hpp>
#include <snap/meshblock_impl.hpp>

// application
#include <application/application.hpp>

// harp
#include "absorber.hpp"

Absorber::Absorber(std::string name) : name_(name) {
  Application::Logger app("harp");
  app->Log("Create Absorber " + name_);
}

Absorber::Absorber(MeshBlock* pmb, std::string name,
                   std::vector<std::string> species,
                   std::map<std::string, Real> params)
    : name_(name), params_(params) {
  Application::Logger app("harp");
  app->Log("Create Absorber " + name_);

  for (auto s : species) {
    imols_.push_back(pmb->pimpl->GetSpeciesId(s));
  }

  app->Log("Dependent species ids = ", imols_);
}

Absorber::~Absorber() {
  Application::Logger app("harp");
  app->Log("Destroy Absorber " + name_);
}

void Absorber::LoadCoefficient(std::string fname, size_t bid) {}

Real Absorber::GetAttenuation(Real wave1, Real wave2,
                              CellVariables const& var) const {
  return 0.;
}

Real Absorber::GetSingleScatteringAlbedo(Real wave1, Real wave2,
                                         CellVariables const& var) const {
  return 0.;
}

void Absorber::GetPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                CellVariables const& var, int np) const {}
