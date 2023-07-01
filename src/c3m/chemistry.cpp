// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// canoe
#include <configure.hpp>

// chemistry
#include "chemistry.hpp"

Chemistry::Chemistry(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NCHEMISTRY == 0) return;

  r.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD + NTRACER, NCHEMISTRY);
  s.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD + NTRACER, NCHEMISTRY);
}
