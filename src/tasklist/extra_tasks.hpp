#ifndef SRC_TASKLIST_EXTRA_TASKS_HPP_
#define SRC_TASKLIST_EXTRA_TASKS_HPP_

// athena
#include <athena/globals.hpp>

// forward declarations
class Mesh;
class MeshBlock;
class TaskID;
class TaskList;
class ParameterInput;

class InversionTasks : public TaskList {
 public:
  InversionTaskList(ParameterInput *pin, Mesh *pm);
  ~InversionTaskList() {}
  // void AddInversionTask(uint64_t id, uint64_t dep);

  TaskStatus CalculateGradient(MeshBlock *pmb, int step);
  TaskStatus Sample(MeshBlock *pmb, int step);
  TaskStatus Optimize(MeshBlock *pmb, int step);

 private:
  void AddTask(const TaskID &id, const TaskID &dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

namespace InversionTaskNames {
const TaskID NONE(0);
const TaskID CALC_GRAD(1);
const TaskID OPTIMIZE(2);
const TaskID SAMPLE(3);
}  // namespace InversionTaskNames

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

#endif  // SRC_TASKLIST_EXTRA_TASKS_HPP_
