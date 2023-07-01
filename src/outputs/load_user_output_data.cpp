// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>

// canoe
#include <impl.hpp>

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
      pod->name = pband->GetName() + "tau";
      pod->data.InitWithShallowSlice(pband->btau, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_ += 1;
    }
  }

  if (output_params.variable.compare("rad") == 0 ||
      output_params.variable.compare("radflux") == 0) {
    // flux up and down
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "flxup";
    pod->data.InitWithShallowSlice(prad->flxup, 4, 0, 1);
    AppendOutputDataNode(pod);
    num_vars_ += 1;

    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "flxdn";
    pod->data.InitWithShallowSlice(prad->flxdn, 4, 0, 1);
    AppendOutputDataNode(pod);
    num_vars_ += 1;
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

  // mcmc inversion
  if (output_params.variable.compare("mcmc")) {
    pod = new OutputData;
    AppendOutputDataNode(pod);
    num_vars_ += 1;
  }
}
