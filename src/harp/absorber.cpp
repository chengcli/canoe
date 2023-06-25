// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/cell_variables.hpp>

// application
#include <application/application.hpp>

// harp
#include "absorber.hpp"

Absorber::Absorber(std::string name) : name_(name) {
  Application::Logger app("harp");
  app->Log("Absorber " + name_ + " is created");
}

Absorber::~Absorber() {
  Application::Logger app("harp");
  app->Log("Absorber " + name_ + " is destroyed");
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
