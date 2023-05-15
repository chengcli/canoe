// C/C++ headers
#include <iostream>
#include <stdexcept>
#include <ctime>

// Athena++ headers
#include <globals.hpp>
#include <utils/utils.hpp>
#include <mesh/mesh.hpp>
#include <utils/utils.hpp>

// debugger headers
#include <debugger.hpp>

// cliutils header
#include <configure.hpp>
#include "command_line.hpp"
#include "program_setup.hpp"

#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

namespace Globals {
  int           mpi_tag_ub;
  clock_t       tstart;
  std::uint64_t mbcnt;
  CommandLine   *cli;
}


void program_start(int argc, char **argv) {
  Globals::cli = new CommandLine(argc, argv);

  Globals::tstart = clock();

  SignalHandler::SignalHandlerInit();
  if (Globals::my_rank == 0 && CommandLine::wtlim > 0)
    SignalHandler::SetWallTimeAlarm(CommandLine::wtlim);

#ifdef MPI_PARALLEL
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    std::runtime_error("MPI initialization failed");
  }

  // Get process id (rank) in MPI_COMM_WORLD
  if (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &(Globals::my_rank))) {
    std::runtime_error("MPI_Comm_rank failed");
  }

  // Get total number of MPI processes (ranks)
  if (MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD, &Globals::nranks)) {
    std::runtime_error("MPI_Comm_size failed");
  }

  // Get maximum value of MPI tag
  MPI_Aint* tag_ub_ptr;
  int att_flag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub_ptr, &att_flag);
  Globals::mpi_tag_ub = *tag_ub_ptr;

#else // no MPI
  Globals::my_rank    = 0;
  Globals::nranks     = 1;
  Globals::mpi_tag_ub = 0;
#endif

  if (Globals::my_rank == 0) {
    std::cout << std::endl;
    std::cout << "######################################################" << std::endl;
    std::cout << "##                 EXDACS MODEL                     ##" << std::endl;
    std::cout << "######################################################" << std::endl;
  }

  pdebug = std::make_unique<Debugger>(1);
  Globals::mbcnt = 0;
  pdebug->Enter("Main");
}

void program_end()
{
  clock_t tstop = clock();
  double cpu_time = (tstop>Globals::tstart ? static_cast<double> (tstop-Globals::tstart) :
                     1.0)/static_cast<double> (CLOCKS_PER_SEC);
  std::cout << "cpu time used  = " << cpu_time << std::endl;
}

void program_end(Mesh *pmesh)
{
  delete Globals::cli;

  clock_t tstop = clock();

  if (Globals::my_rank == 0) {
    pmesh->OutputCycleDiagnostics();
    if (SignalHandler::GetSignalFlag(SIGTERM) != 0) {
      std::cout << std::endl << "Terminating on Terminate signal" << std::endl;
    } else if (SignalHandler::GetSignalFlag(SIGINT) != 0) {
      std::cout << std::endl << "Terminating on Interrupt signal" << std::endl;
    } else if (SignalHandler::GetSignalFlag(SIGALRM) != 0) {
      std::cout << std::endl << "Terminating on wall-time limit" << std::endl;
    } else if (pmesh->ncycle == pmesh->nlim) {
      std::cout << std::endl << "Terminating on cycle limit" << std::endl;
    } else {
      std::cout << std::endl << "Terminating on time limit" << std::endl;
    }

    std::cout << "time=" << pmesh->time << " cycle=" << pmesh->ncycle << std::endl;
    std::cout << "tlim=" << pmesh->tlim << " nlim=" << pmesh->nlim << std::endl;

    if (pmesh->adaptive) {
      std::cout << std::endl << "Number of MeshBlocks = " << pmesh->nbtotal
                << "; " << pmesh->nbnew << "  created, " << pmesh->nbdel
                << " destroyed during this simulation." << std::endl;
    }

    // Calculate and print the zone-cycles/cpu-second and wall-second
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;
#endif
    clock_t tstop = clock();
    double cpu_time = (tstop>Globals::tstart ? static_cast<double> (tstop-Globals::tstart) :
                       1.0)/static_cast<double> (CLOCKS_PER_SEC);
    std::uint64_t zonecycles =
        Globals::mbcnt*static_cast<std::uint64_t> (pmesh->my_blocks(0)->GetNumberOfMeshBlockCells());
    double zc_cpus = static_cast<double> (zonecycles) / cpu_time;

    std::cout << std::endl << "zone-cycles = " << zonecycles << std::endl;
    std::cout << "cpu time used  = " << cpu_time << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    double zc_omps = static_cast<double> (zonecycles) / omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
  }

  if (Globals::my_rank == 0 && CommandLine::wtlim > 0)
    SignalHandler::CancelWallTimeAlarm();

  pdebug->Leave();
}
