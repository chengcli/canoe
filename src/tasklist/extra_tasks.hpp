#ifndef SRC_TASKLIST_EXTRA_TASKS_HPP_
#define SRC_TASKLIST_EXTRA_TASKS_HPP_

// athena
#include <athena/globals.hpp>
#include <athena/task_list/task_list.hpp>

// forward declarations
class Mesh;
class MeshBlock;
class TaskID;
class TaskList;
class ParameterInput;

class InversionTasks : public TaskList {
 public:
  InversionTasks(ParameterInput *pin, Mesh *pm);
  ~InversionTasks() {}
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
  TaskStatus AddSourceTerms(MeshBlock *pmb, int stage);

  TaskStatus AddFluxToConserved(MeshBlock *pmb, int stage);
  TaskStatus ImplicitCorrection(MeshBlock *pmb, int stage);
  TaskStatus UpdateAllConserved(MeshBlock *pmb, int stage);

  // turbulence tasks
  TaskStatus CalculateTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus SendTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus ReceiveTurbulenceFlux(MeshBlock *pmb, int stage);
  TaskStatus IntegrateTurbulence(MeshBlock *pmb, int stage);
  TaskStatus SendTurbulence(MeshBlock *pmb, int stage);
  TaskStatus ReceiveTurbulence(MeshBlock *pmb, int stage);
  TaskStatus SetBoundariesTurbulence(MeshBlock *pmb, int stage);

 protected:
  void AddTask(TaskID const &id, TaskID const &dep) override;
};

//! This should track the largest task ID in the athena/task_list/task_list.hpp
namespace HydroIntegratorTaskNames {
// implicit
const TaskID IMPLICIT_CORR(80);
const TaskID ADD_FLX_CONS(81);
const TaskID UPDATE_ALLCONS(82);

// turbulence
const TaskID CALC_TURBFLX(83);
const TaskID SEND_TURBFLX(84);
const TaskID RECV_TURBFLX(85);
const TaskID INT_TURB(86);
const TaskID SEND_TURB(87);
const TaskID RECV_TURB(88);
const TaskID SETB_TURB(89);
}  // namespace HydroIntegratorTaskNames

#endif  // SRC_TASKLIST_EXTRA_TASKS_HPP_
