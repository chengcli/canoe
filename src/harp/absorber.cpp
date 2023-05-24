#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <debugger/debugger.hpp>
#include <snap/cell_variables.hpp>
// #include <snap/mesh/meshblock_impl.hpp>
// #include "../radiation/radiation.hpp"
#include "absorber.hpp"

extern std::unique_ptr<Debugger> pdebug;

Absorber::Absorber(MeshBlock* pmb, ParameterInput* pin, std::string bname,
                   std::string name)
    : name_(name) {
  pdebug->Enter("Absorber " + name);
  pdebug->Leave();
}

Absorber::~Absorber() {}

void Absorber::loadCoefficient(std::string fname, int bid) {}

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
