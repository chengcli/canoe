// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// canoe
#include <configure.hpp>

// dusts
#include "cloud.hpp"

Cloud::Cloud(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NCLOUD == 0) return;

  w.InitWithShallowSlice(pmb->pscalars->r, 4, 0, NCLOUD);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, 0, NCLOUD);
}
