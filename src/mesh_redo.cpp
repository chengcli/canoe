// athena
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

void Mesh::SaveAllStates()
{
  for (int b = 0; b < nblocal; ++b) {
    my_blocks(b)->pimpl->SaveAllStates();
  }
}

void Mesh::LoadAllStates()
{
  for (int b = 0; b < nblocal; ++b) {
    my_blocks(b)->pimpl->LoadAllStates();
  }
}

void Mesh::DecreaseTimeStep()
{
  dt /= 2;
}

bool Mesh::CheckAllValid() const
{
  bool valid = true;
  for (int b = 0; b < nblocal; ++b)
    if (!my_blocks(b)->pimpl->IsStateValid())
      valid = false;

  bool all_valid = valid;

#ifdef MPI_PARALLEL
  // MPI call to check if all blocks are valid
  MPI_Allreduce(&valid, &all_valid, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
#endif  // MPI_PARALLEL
  
  return all_valid;
}
