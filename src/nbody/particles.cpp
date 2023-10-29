//! \file particles.cpp
//! \brief Implements class ParticleBase and helper functions

// C++ headers
#include <string>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// Athena++
#include <athena/mesh/mesh.hpp>

// utils
#include <utils/vectorize.hpp>

// n-body
#include "particle_data.hpp"
#include "particles.hpp"

// constructor, initializes data structure and parameters
ParticleBase::ParticleBase(MeshBlock *pmb, std::string name)
    : NamedGroup(name), linked_flag_(false), pmy_block_(pmb) {
  Application::Logger app("n-body");
  app->Log("Initialize ParticleBase");

  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  weight.NewAthenaArray(nc3, nc2, nc1);
  charge.NewAthenaArray(nc3, nc2, nc1);
  pd_in_cell_.NewAthenaArray(nc3, nc2, nc1);
}

ParticleBase::~ParticleBase() {
  Application::Logger app("n-body");
  app->Log("Destroy ParticleBase");
}

AllParticles ParticlesFactory::CreateAllParticles(MeshBlock *pmb,
                                                  ParameterInput *pin) {
  AllParticles all_particles;

  std::string str = pin->GetOrAddString("particles", "particles", "");
  auto particle_names = Vectorize<std::string>(str.c_str());

  for (auto const &name : particle_names) {
    if (name == "2pcp") {
      all_particles.push_back(std::make_shared<ParticleBase>(pmb, name));
    } else if (name == "scp") {
      all_particles.push_back(std::make_shared<ParticleBase>(pmb, name));
    } else {
      throw NotImplementedError("Particle '" + name + "' unrecognized.");
    }
  }

  return all_particles;
}

namespace ParticlesHelper {

ParticlePtr find_particle(AllParticles const &pts, std::string name) {
  for (auto &pt : pts) {
    if (pt->GetName() == name) return pt;
  }
  return nullptr;
}

} // namespace ParticlesHelper

