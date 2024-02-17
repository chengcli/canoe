// C/C++
#include <iostream>
#include <memory>

// canoe
#include <configure.hpp>

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

// athena
#include <athena/parameter_input.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/task_list/task_list.hpp>

// tasklist
#include <tasklist/task_list_factory.hpp>

// application
#include <application/command_line.hpp>
#include <application/signal.hpp>
#include <application/application.hpp>
#include <application/exceptions.hpp>

namespace Globals {

const char* search_paths = "@CMAKE_SOURCE_DIR@/src/snap/thermodynamics/chemdata/:@CMAKE_SOURCE_DIR@/data/";

const char* banner =
    "\
######################################################\n\
##              CANOE MODEL STARTS                  ##\n\
######################################################";

}  // namespace Globals

void mesh_setup(ParameterInput* &pin, Mesh* &mesh);
void mesh_destroy(ParameterInput* &pin, Mesh* &mesh, int mbcnt);

using OutputsPtr = std::unique_ptr<Outputs>;

int main(int argc, char **argv)  {
  Application::Start(argc, argv);

#ifdef ENABLE_GLOG
  Application::Logger app("main");

  std::string prog_name = argv[0];

  // Extract the base name
  size_t last_slash_pos = prog_name.find_last_of("/");
  std::string base_name = (last_slash_pos == std::string::npos)
                              ? prog_name
                              : prog_name.substr(last_slash_pos + 1);
  std::string log_dir_name = base_name + ".glog";

  struct stat st = {0};
  if (stat(log_dir_name.c_str(), &st) == -1) { // Check if the directory doesn't exist.
      if (mkdir(log_dir_name.c_str(), 0755) == -1) { // Mode 0755 gives rwx permissions for everyone.
          perror("Error creating directory");
      } else {
          app->Log(log_dir_name + " created.");
      }
  } else {
    app->Log(log_dir_name + " exists.");
  }

  google::InitGoogleLogging(argv[0]);
#endif

  auto cli = CommandLine::GetInstance();
  auto sig = Signal::GetInstance();

  ParameterInput *pinput;
  Mesh *pmesh;

  mesh_setup(pinput, pmesh);

  // CONSTRUCT AND INITIALIZE TASKlIST
  auto ptlist = TaskListFactory::CreateFrom(pinput, pmesh);

#ifdef ENABLE_GLOG
  FLAGS_log_dir = log_dir_name;
  char log_name[FILENAME_MAX];
  snprintf(log_name, FILENAME_MAX, "%s.%d", base_name.c_str(), Globals::my_rank);
  google::SetLogSymlink(google::GLOG_INFO, log_name);
  google::SetLogSymlink(google::GLOG_WARNING, log_name);
  google::SetLogSymlink(google::GLOG_ERROR, log_name);
#endif

  // CHANGE TO RUN DIRECTORY, INITIALIZE OUTPUTS OBJECT, AND MAKE OUTPUT OF ICs
  OutputsPtr pouts;
  try {
    //ChangeRunDir(cli->prundir);
    pouts = std::make_unique<Outputs>(pmesh, pinput);
    if (cli->res_flag == 0)
      pouts->MakeOutputs(pmesh, pinput);
  }
  catch(std::bad_alloc& ba) {
    throw RuntimeError("main", "Output memory allocation failed");
  }
  catch(std::exception const& ex) {
    throw RuntimeError("main", ex.what());
  }

  // START OF MAIN INTEGRATION LOOP
  // For performance, there is no error handler protecting this step (except outputs)

  if (Globals::my_rank == 0) {
    std::cout << "\nSetup complete, entering main loop...\n";
  }

  std::uint64_t mbcnt = 0;
  int attempts = 0;

  while ((pmesh->time < pmesh->tlim) &&
         (pmesh->nlim < 0 || pmesh->ncycle < pmesh->nlim)) {
    if (Globals::my_rank == 0)
      pmesh->OutputCycleDiagnostics();

    for (int stage=1; stage<=ptlist->nstages; ++stage) {
      ptlist->DoTaskListOneStage(pmesh, stage);
    }

    pmesh->UserWorkInLoop();

    pmesh->ncycle++;
    pmesh->time += pmesh->dt;
    mbcnt += pmesh->nbtotal;
    pmesh->step_since_lb++;

    pmesh->LoadBalancingAndAdaptiveMeshRefinement(pinput);

    pmesh->NewTimeStep();

    try {
      if (pmesh->time < pmesh->tlim) // skip the final output as it happens later
        pouts->MakeOutputs(pmesh,pinput);
    }
    catch(std::bad_alloc& ba) {
      throw RuntimeError("main", "MakeOutput memory allocation failed");
    }
    catch(std::exception const& ex) {
#ifdef MPI_PARALLEL
      MPI_Finalize();
#endif
      throw RuntimeError("main", ex.what());
    }

    // check for signals
    if (sig->CheckSignalFlags() != 0) {
      break;
    }
  }

  // END OF MAIN INTEGRATION LOOP
  // Make final outputs, print diagnostics, clean up and terminate

  // Output the final cycle diagnostics and make the final outputs
  if (Globals::my_rank == 0)
    pmesh->OutputCycleDiagnostics();

  pmesh->UserWorkAfterLoop(pinput);

  try {
    pouts->MakeOutputs(pmesh,pinput,true);
  }
  catch(std::bad_alloc& ba) {
    throw RuntimeError("main", "MakeOutput memory allocation failed");
  }
  catch(std::exception const& ex) {
#ifdef MPI_PARALLEL
    MPI_Finalize();
#endif
    throw RuntimeError("main", ex.what());
  }

  // Print diagnostic messages related to the end of the simulation
  mesh_destroy(pinput, pmesh, mbcnt);

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif

  Application::Destroy();
}
