// athena
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/io_wrapper.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/command_line.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.h>

#include "impl.hpp"
#include "index_map.hpp"

// snap
#include "snap/thermodynamics/thermodynamics.hpp"

// n-body
// #include "nbody/particle_data.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

void mesh_setup(ParameterInput*& pinput, Mesh*& pmesh) {
  IOWrapper infile, restartfile;
  auto cli = CommandLine::GetInstance();

  try {
    pinput = new ParameterInput;

    if (cli->res_flag == 1) {
      restartfile.Open(cli->restart_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(restartfile);
      // If both -r and -i are specified, make sure next_time gets corrected.
      // This needs to be corrected on the restart file because we need the old
      // dt.
      if (cli->iarg_flag == 1) pinput->RollbackNextTime();
      // leave the restart file open for later use
    }

    if (cli->iarg_flag == 1) {
      // if both -r and -i are specified, override the parameters using the
      // input file
      infile.Open(cli->input_filename, IOWrapper::FileMode::read);
      pinput->LoadFromFile(infile);
      infile.Close();
    }
    pinput->ModifyFromCmdline(cli->argc, cli->argv);
  } catch (std::bad_alloc& ba) {
    if (cli->res_flag == 1) restartfile.Close();
    throw RuntimeError(
        "main", "memory allocation failed initializing class ParameterInput:");
  } catch (std::exception const& ex) {
    if (cli->res_flag == 1) restartfile.Close();
    throw RuntimeError("main", ex.what());
  }

  // index map
  IndexMap::InitFromAthenaInput(pinput);

  // thermodynamics
  Thermodynamics::InitFromAthenaInput(pinput);

  // n-body
  // ParticleHelper::commit_mpi_particle_data();

  try {
    if (cli->res_flag == 0) {
      pmesh = new Mesh(pinput, cli->mesh_flag);
    } else {
      pmesh = new Mesh(pinput, restartfile, cli->mesh_flag);
    }
  } catch (std::bad_alloc& ba) {
    if (cli->res_flag == 1) restartfile.Close();
    throw RuntimeError("main",
                       "memory allocation failed initializing class Mesh:");
  } catch (std::exception const& ex) {
    if (cli->res_flag == 1) restartfile.Close();
    throw RuntimeError("main", ex.what());
  }

  // With current mesh time possibly read from restart file, correct next_time
  // for outputs
  if (cli->iarg_flag == 1 && cli->res_flag == 1) {
    // if both -r and -i are specified, ensure that next_time  >= mesh_time - dt
    pinput->ForwardNextTime(pmesh->time);
  }

  // Dump input parameters and quit if code was run with -n option.
  if (cli->narg_flag) {
    if (Globals::my_rank == 0) pinput->ParameterDump(std::cout);
    if (cli->res_flag == 1) restartfile.Close();
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    exit(0);
  }

  // everything works fine here, close the restart file
  if (cli->res_flag == 1) restartfile.Close();

  // Quit if -m was on cmdline.  This option builds and outputs mesh structure.
  if (cli->mesh_flag > 0) {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    exit(0);
  }

  // set up additional components
  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock* pmb = pmesh->my_blocks(b);
    pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
  }

  // initialize mesh
  pmesh->Initialize(cli->res_flag, pinput);
}
