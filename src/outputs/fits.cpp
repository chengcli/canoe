/*! \file fits.cpp
 *  \brief writes output data in FITS format (for inversion tasks)
 */

// C/C++
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/user_outputs.hpp>

// canoe
#include <impl.hpp>

// inversion
#include <inversion/inversion.hpp>

#ifdef FITSOUTPUT
extern "C" {
#include <fitsio.h>
}

FITSOutput::FITSOutput(OutputParameters oparams) : OutputType(oparams) {}

void FITSOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  // TODO : saftgard for pfit == null
  // create filename: "file_basename"+"."+"file_id"+"."+XXXXX+".fits",
  // where XXXXX = 5-digit file_number
  if (output_params.file_number > 0) {
    std::string fname = "!";  // clobber
    char number[6];
    int err;
    snprintf(number, sizeof(number), "%05d", output_params.file_number);

    fname.append(output_params.file_basename);
    fname.append(".");
    fname.append(output_params.file_id);
    fname.append(".");
    fname.append(number);
    fname.append(".fits");

    /*for (int i = 0; i < pm->nblocal; ++i) {
      MeshBlock *pmb = pm->my_blocks(i);
      auto &pfit = pmb->pimpl->fitq.back();
      pfit->MakeMCMCOutputs(fname);
      pfit->ResetChain();
    }*/
  }

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number",
                  output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
}

#endif  // FITSOUTPUT
