#ifndef IMPLICIT_HYDRO_TASKS_HPP
#define IMPLICIT_HYDRO_TASKS_HPP

// canoe headers
#include <configure.hpp>

// athena
#include <task_list/task_list.hpp>

// snap
#include "../constants.hpp"

class ParameterInput;
class Mesh;

class ImplicitHydroTasks : public TimeIntegratorTaskList {
 public:
  ImplicitHydroTasks(ParameterInput *pin, Mesh *pm);

  // implicit tasks
  TaskStatus IntegrateHydro(MeshBlock *pmb, int stage);
  TaskStatus UpdateHydro(MeshBlock *pmb, int stage);
  TaskStatus AddSourceTerms(MeshBlock *pmb, int stage);

 protected:
  void AddTask(TaskID const &id, TaskID const &dep) override;
};

namespace HydroIntegratorTaskNames {
const TaskID UPDATE_HYD(Globals::next_task_id);
}

#endif
