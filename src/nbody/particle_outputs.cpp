//! \file particle_outputs.cpp
//! \brief Output functions for particles

// athena
#include <athena/outputs/outputs.hpp>

// n-body
#include "particles.hpp"

bool ParticleBase::ShouldMeshOutput(std::string variable_name) const {
  return variable_name.compare("particle") == 0;
}

void ParticleBase::LoadMeshOutputData(OutputType *out, int *num_vars) const {
  OutputData *pod = new OutputData;
  pod->type = "SCALARS";
  pod->name = "p-weight";
  pod->data.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(weight), 4, 0,
                                 1);
  out->AppendOutputDataNode(pod);
  *num_vars += 1;

  pod = new OutputData;
  pod->type = "SCALARS";
  pod->name = "p-charge";
  pod->data.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(charge), 4, 0,
                                 1);
  out->AppendOutputDataNode(pod);
  *num_vars += 1;
}
