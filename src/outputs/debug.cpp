//! \file debug.cpp
//  \brief writes debug output data, max and min quantities and their locations that are output
//         frequently in time to trace extreme values.

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cfloat>

#include <athena/athena.hpp>
#include <athena/globals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe headers
#include <configure.hpp>
#include "user_outputs.hpp"

//----------------------------------------------------------------------------------------
//! \fn void DebugOutput::WriteOutputFile()
//  \brief Writes a history file

void DebugOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  struct {
    Real value;
    int rank;
  } vmin[NHYDRO], vmax[NHYDRO];

  struct {
    Real x1, x2, x3;
  } xmin[NHYDRO], xmax[NHYDRO];

  for (int n = 0; n < NHYDRO; ++n) {
    vmin[n].value = FLT_MAX;
    vmin[n].rank = Globals::my_rank;
    vmax[n].value = FLT_MIN;
    vmax[n].rank = Globals::my_rank;
    xmin[n].x1 = xmin[n].x2 = xmin[n].x3 = FLT_MAX;
    xmax[n].x1 = xmax[n].x2 = xmax[n].x3 = FLT_MIN;
  }

  // Loop over MeshBlocks
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    Hydro *phydro = pmb->phydro;
    Coordinates *pcoord = pmb->pcoord;

    int is = pmb->is, js = pmb->js, ks = pmb->ks;
    int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
    //out_is=pmb->is; out_ie=pmb->ie;
    //out_js=pmb->js; out_je=pmb->je;
    //out_ks=pmb->ks; out_ke=pmb->ke;
    // ghost cells are included in the calculations
    //out_is -= NGHOST; out_ie += NGHOST;
    //if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
    //if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}

    // calculate maximum and minimum values over cells
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j= js; j<= je; ++j)
          for (int i= is; i<= ie; ++i) {
            if (phydro->w(n,k,j,i) < vmin[n].value) {
              vmin[n].value = phydro->w(n,k,j,i);
              xmin[n].x1 = pcoord->x1v(i);
              xmin[n].x2 = pcoord->x2v(j);
              xmin[n].x3 = pcoord->x3v(k);
            }
            if (phydro->w(n,k,j,i) > vmax[n].value) {
              vmax[n].value = phydro->w(n,k,j,i);
              xmax[n].x1 = pcoord->x1v(i);
              xmax[n].x2 = pcoord->x2v(j);
              xmax[n].x3 = pcoord->x3v(k);
            }
          }
  }

#ifdef MPI_PARALLEL
  // gather all nodes and synchronize
  // TODO: xmin, xmax semm not correct
  MPI_Allreduce(MPI_IN_PLACE, vmin, NHYDRO, MPI_REAL_INT, MPI_MINLOC, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, vmax, NHYDRO, MPI_REAL_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  // distribute the coordinates of the extreme values to master node
  MPI_Status status;
  MPI_Request request[2*NHYDRO];
  for (int n = 0; n < 2*NHYDRO; ++n)
    request[n] = MPI_REQUEST_NULL;

  // send and receive coordinate for min values
  for (int n = 0; n < NHYDRO; ++n) {
    if (Globals::my_rank == 0) {
      if (Globals::my_rank == vmin[n].rank)
        continue;
      else
        MPI_Irecv(xmin + n, 3, MPI_ATHENA_REAL, vmin[n].rank, n, MPI_COMM_WORLD, request+n);
    } else {
      if (Globals::my_rank == vmin[n].rank)
        MPI_Isend(xmin + n, 3, MPI_ATHENA_REAL, 0, n, MPI_COMM_WORLD, request+n);
      else
        continue;
    }
  }

  // send and receive coordinate for max values
  for (int n = 0; n < NHYDRO; ++n) {
    if (Globals::my_rank == 0) {
      if (Globals::my_rank == vmax[n].rank)
        continue;
      else
        MPI_Irecv(xmax + n, 3, MPI_ATHENA_REAL, vmax[n].rank, n+NHYDRO, MPI_COMM_WORLD, request+NHYDRO+n);
    } else {
      if (Globals::my_rank == vmax[n].rank)
        MPI_Isend(xmax + n, 3, MPI_ATHENA_REAL, 0, n+NHYDRO, MPI_COMM_WORLD, request+NHYDRO+n);
      else
        continue;
    }
  }

  // blocks and waits for master node to receive all data
  for (int n = 0; n < 2*NHYDRO; ++n)
    MPI_Wait(request+n, &status);
#endif

  // only the master rank writes the file
  // create filename: "file_basename" + ".hst".  There is no file number.
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign(output_params.file_basename);
    fname.append(".dbg");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen(fname.c_str(),"a")) == NULL){
      msg << "### FATAL ERROR in function [OutputType::DebugFile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      ATHENA_ERROR(msg);
    }

    // If this is the first output, write header
    int iout = 1;
    if (output_params.file_number == 0) {
      fprintf(pfile,"# Athena++ debug data\n"); // descriptor is first line
      fprintf(pfile,"# [%d]=time    ", iout++);
      fprintf(pfile,"[%d]=id       ", iout++);
      fprintf(pfile,"[%d]=dmin     ", iout++);
      fprintf(pfile,"[%d]=dmax     ", iout++);
      fprintf(pfile,"[%d]=x1min    ", iout++);
      fprintf(pfile,"[%d]=x1max    ", iout++);
      fprintf(pfile,"[%d]=x2min    ", iout++);
      fprintf(pfile,"[%d]=x2max    ", iout++);
      fprintf(pfile,"[%d]=x3min    ", iout++);
      fprintf(pfile,"[%d]=x3max    ", iout++);
      fprintf(pfile,"\n"); // terminate line
    }

    // write debug variables
    for (int n = 0; n < NHYDRO; ++n) {
      if (n == 0)
        fprintf(pfile, output_params.data_format.c_str(), pm->time);
      else
        fprintf(pfile,"             ");
      fprintf(pfile, " --  %2d  -- ", n);
      fprintf(pfile, output_params.data_format.c_str(), vmin[n].value);
      fprintf(pfile, output_params.data_format.c_str(), vmax[n].value);
      fprintf(pfile, output_params.data_format.c_str(), xmin[n].x1);
      fprintf(pfile, output_params.data_format.c_str(), xmax[n].x1);
      fprintf(pfile, output_params.data_format.c_str(), xmin[n].x2);
      fprintf(pfile, output_params.data_format.c_str(), xmax[n].x2);
      fprintf(pfile, output_params.data_format.c_str(), xmin[n].x3);
      fprintf(pfile, output_params.data_format.c_str(), xmax[n].x3);
      fprintf(pfile,"\n");
    }
    fclose(pfile);
  }

  // increment counters, clean up
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
