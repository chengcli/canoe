#ifndef NBODY_PARTICLE_EXCHANGER_HPP_
#define NBODY_PARTICLE_EXCHANGER_HPP_

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// C/C++ header
#include <vector>

// communication
#include <communicator/exchanger.hpp>

// n-body
#include "particle_data.hpp"
#include "particle_group.hpp"

//! Defines the class for managing buffers for transporting particles.
class ParticleExchanger : public Exchanger<ParticleData> {
 protected:
  Particle const *pmy_particle;

 public:
  // functions
  ParticleExchanger(Particle *pg) : pmy_particle(pg) {}

  ~ParticleExchanger() {}

  void DetachTo(std::vector<ParticleData> &buffer) override;
  bool AttachTo(std::vector<ParticleData> &container) const override;

 protected:
  bool checkOutOfDomain(ParticleData &p) const;

  MeshBlock const *getMeshBlock() const override {
    return pmy_particle->pmy_block;
  }
};

#endif  // NBODY_PARTICLE_EXCHANGER_HPP_
