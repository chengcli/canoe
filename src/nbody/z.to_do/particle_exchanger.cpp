// C/C++ headers
#include <sstream>
// #include <cstddef>
#include <functional>  // hash
#include <iostream>

// n-body
#include "particle_data.hpp"
#include "particle_exchanger.hpp"
#include "particle_group.hpp"

#ifdef MPI_PARALLEL
// defined in main.cpp
extern MPI_Datatype MPI_PARTICLE;
#endif

ParticleExchanger::ParticleExchanger(ParticleGroup *pg) : pmy_particle(pg) {}

ParticleExchanger::~ParticleExchanger() {}

void ParticleExchanger::SendParticle() {}

void ParticleExchanger::RecvParticle() {}

void ParticleExchanger::ClearBoundary() {}
