// C/C++
#include <iostream>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/io_wrapper.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/globals.hpp>

// canoe
#include <configure.h>

#include "index_map.hpp"

// outputs
#include "outputs/output_utils.hpp"

void mesh_destroy(ParameterInput *&pinput, Mesh *&pmesh, int mbcnt) {
  clock_t tstop = clock();

  if (Globals::my_rank == 0) {
    pmesh->OutputCycleDiagnostics();
    if (pmesh->ncycle == pmesh->nlim) {
      std::cout << std::endl << "Terminating on cycle limit" << std::endl;
    }

    std::cout << "time=" << pmesh->time << " cycle=" << pmesh->ncycle
              << std::endl;
    std::cout << "tlim=" << pmesh->tlim << " nlim=" << pmesh->nlim << std::endl;

    if (pmesh->adaptive) {
      std::cout << std::endl
                << "Number of MeshBlocks = " << pmesh->nbtotal << "; "
                << pmesh->nbnew << "  created, " << pmesh->nbdel
                << " destroyed during this simulation." << std::endl;
    }

    // Calculate and print the zone-cycles/cpu-second and wall-second
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;
#endif
    clock_t tstop = clock();
    double cpu_time =
        (tstop > Globals::tstart ? static_cast<double>(tstop - Globals::tstart)
                                 : 1.0) /
        static_cast<double>(CLOCKS_PER_SEC);
    std::uint64_t zonecycles =
        mbcnt * static_cast<std::uint64_t>(
                    pmesh->my_blocks(0)->GetNumberOfMeshBlockCells());
    double zc_cpus = static_cast<double>(zonecycles) / cpu_time;

    std::cout << std::endl << "zone-cycles = " << zonecycles << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    double zc_omps = static_cast<double>(zonecycles) / omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
  }

  delete pinput;
  delete pmesh;

  IndexMap::Destroy();
  MetadataTable::Destroy();
}
