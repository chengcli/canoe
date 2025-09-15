// athena
#include <athena/athena.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/user_outputs.hpp>
#include <athena/scalars/scalars.hpp>

// kintera
#include <kintera/utils/serialize.hpp>

// canoe
#include <interface/hydro.hpp>

#include "output_utils.hpp"

PTOutput::PTOutput(OutputParameters oparams) : OutputType(oparams) {}

void PTOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // create filename: "file_basename"+"."+"file_id"+"."+XXXXX+".pt",
  // where XXXXX = 5-digit file_number
  auto pmeta = MetadataTable::GetInstance();

  std::string fname;
  char number[6];
  int err;
  snprintf(number, sizeof(number), "%05d", output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".");
  fname.append(std::to_string(Globals::my_rank));
  fname.append(".pt");

  // Loop over MeshBlocks
  std::map<std::string, torch::Tensor> data;
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);

    data["hydro_u/" + std::to_string(b)] = get_all(pmb->phydro->u);
    data["scalar_s/" + std::to_string(b)] = get_all(pmb->pscalars->s);
  }  // end loop over MeshBlocks

  kintera::save_tensors(data, fname);

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number",
                  output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
}
