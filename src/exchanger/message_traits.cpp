// C/C++
#include <functional>
#include <string>

// athena++
#include <athena/athena.hpp>
#include <athena/bvals/bvals.hpp>

// nbody
#include <nbody/particles.hpp>

// harp
#include <harp/radiation_band.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>

MPI_Datatype MPI_PARTICLE_DATA;
MPI_Datatype MessageTraits<ParticleBase>::mpi_type = MPI_PARTICLE_DATA;
MPI_Datatype MessageTraits<RadiationBand>::mpi_type = MPI_ATHENA_REAL;

#endif  // MPI_PARALLEL

namespace MessageHelper {
int mpi_tag_ub;

int create_mpi_tag(int lid, int tid, std::string name) {
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, 0);

  std::string str = name + std::to_string(tag);
  return std::hash<std::string>{}(str) % (mpi_tag_ub);
}
}  // namespace MessageHelper
