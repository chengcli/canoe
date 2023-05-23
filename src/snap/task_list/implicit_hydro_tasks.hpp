#ifndef IMPLICIT_HYDRO_TASKS_HPP
#define IMPLICIT_HYDRO_TASKS_HPP

// Athena++ headers
#include <task_list/task_list.hpp>

// canoe headers
#include <configure.hpp>

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
