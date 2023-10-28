// C/C++
#include <iostream>
#include <sstream>
#include <stdexcept>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <impl.hpp>

// inversion
#include <inversion/inversion.hpp>

// tasklist
#include "extra_tasks.hpp"

using namespace InversionTaskNames;

InversionTasks::InversionTasks(ParameterInput *pin, Mesh *pm) {
  nstages = 1;
  std::string task = pin->GetString("inversion", "task");

  // Now assemble list of tasks for each step of inversion task

  if (task == "atm_profile") {
    AddTask(SAMPLE, NONE);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in InversionTasks::InversionTasks" << std::endl
        << "Unrecognized inversion task";
    ATHENA_ERROR(msg);
    // AddTask(CALC_GRAD,NONE);
    // AddTask(OPTIMIZE,CALC_GRAD);
  }
}

void InversionTasks::AddTask(const TaskID &id, const TaskID &dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  if (id == CALC_GRAD) {
    task_list_[ntasks].TaskFunc =
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock *, int)>(
            &InversionTasks::CalculateGradient);
  } else if (id == OPTIMIZE) {
    task_list_[ntasks].TaskFunc =
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock *, int)>(
            &InversionTasks::Optimize);
  } else if (id == SAMPLE) {
    task_list_[ntasks].TaskFunc =
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock *, int)>(
            &InversionTasks::Sample);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in InversionTasks::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void InversionTasks::StartupTaskList(MeshBlock *pmb, int stage) {}

TaskStatus InversionTasks::CalculateGradient(MeshBlock *pmb, int step) {
  // std::cout << "Calculate gradient" << std::endl;
  return TaskStatus::success;
}

TaskStatus InversionTasks::Optimize(MeshBlock *pmb, int step) {
  // std::cout << "Optimize" << std::endl;
  return TaskStatus::success;
}

TaskStatus InversionTasks::Sample(MeshBlock *pmb, int step) {
  // pmb->pimpl->fitq.front()->MCMCMove(pmb->pimpl->prad.get(), pmb->phydro);
  //  std::cout << "Sample" << std::endl;
  return TaskStatus::success;
}
