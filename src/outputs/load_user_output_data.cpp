// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>

// canoe
#include <impl.hpp>
#include <virtual_groups.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// outputs
#include "output_utils.hpp"

void OutputType::loadUserOutputData(MeshBlock *pmb) {
  OutputData *pod;
  auto phyd = pmb->phydro;
  auto prad = pmb->pimpl->prad;

  // vapor
  if (NVAPOR > 0) {
    if (output_params.variable.compare("prim") == 0 ||
        output_params.variable.compare("vapor") == 0) {
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "vapor";
      pod->data.InitWithShallowSlice(phyd->w, 4, 1, NVAPOR);
      AppendOutputDataNode(pod);
      num_vars_ += NVAPOR;
    }

    if (output_params.variable.compare("cons") == 0) {
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "vapor";
      pod->data.InitWithShallowSlice(phyd->u, 4, 1, NVAPOR);
      AppendOutputDataNode(pod);
      num_vars_ += NVAPOR;
    }
  }

  // radiation
  if (output_params.variable.compare("rad") == 0 ||
      output_params.variable.compare("radtau") == 0) {
    for (int b = 0; b < prad->GetNumBands(); ++b) {
      auto pband = prad->GetBand(b);
      // tau
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = pband->GetName() + "-tau";
      pod->data.InitWithShallowSlice(pband->btau, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_ += 1;
    }
  }

  if (output_params.variable.compare("rad") == 0 ||
      output_params.variable.compare("radtime") == 0) {
    // rad time
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "radtime";
    pod->data.InitWithShallowSlice(prad->rtime, 4, 0, 1);
    AppendOutputDataNode(pod);
    num_vars_ += 1;
  }

  if (output_params.variable.compare("rad") == 0 ||
      output_params.variable.compare("radflux") == 0) {
    for (int b = 0; b < prad->GetNumBands(); ++b) {
      auto pband = prad->GetBand(b);
      // flux up
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = pband->GetName() + "-flxup";
      pod->data.InitWithShallowSlice(pband->bflxup, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_ += 1;
    }

    for (int b = 0; b < prad->GetNumBands(); ++b) {
      auto pband = prad->GetBand(b);
      // flux down
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = pband->GetName() + "-flxdn";
      pod->data.InitWithShallowSlice(pband->bflxdn, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_ += 1;
    }
  }

  if (output_params.variable.compare("rad") == 0 ||
      output_params.variable.compare("radtoa") == 0) {
    if (prad->radiance.GetDim3() > 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "radiance";
      pod->data.InitWithShallowSlice(prad->radiance, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_ += 1;

      // for (auto p : prad->bands) p->writeBinRadiance(&output_params);
    }
  }

  /* mcmc inversion
  if (output_params.variable.compare("mcmc")) {
    pod = new OutputData;
    AppendOutputDataNode(pod);
    num_vars_ += 1;
  }*/

  // fits output groups
  for (auto &fits_out : pmb->pimpl->GetFITSOutputGroups()) {
    if (fits_out.lock()->ShouldFITSOutput(output_params.variable))
      fits_out.lock()->LoadFITSOutputData(this, &num_vars_);
  }

  // mesh output groups
  for (auto &mesh_out : pmb->pimpl->GetMeshOutputGroups()) {
    if (mesh_out.lock()->ShouldMeshOutput(output_params.variable))
      mesh_out.lock()->LoadMeshOutputData(this, &num_vars_);
  }
}
