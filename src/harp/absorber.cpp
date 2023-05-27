// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// debugger
#include <debugger/debugger.hpp>

// snap
#include <snap/cell_variables.hpp>

// harp
#include "absorber.hpp"

Absorber::Absorber(MeshBlock* pmb, ParameterInput* pin, std::string bname,
                   std::string name)
    : name_(name) {
  pdebug->Enter("Absorber " + name);
  pdebug->Leave();
}

Absorber::~Absorber() {}

void Absorber::loadCoefficient(std::string fname, size_t bid) {}

Real Absorber::getAttenuation(Real wave1, Real wave2,
                              CellVariables const& var) const {
  return 0.;
}

Real Absorber::getSingleScatteringAlbedo(Real wave1, Real wave2,
                                         CellVariables const& var) const {
  return 0.;
}

void Absorber::getPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                CellVariables const& var, int np) const {}
