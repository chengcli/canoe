#ifndef SRC_SCHEDULERS_HPP_
#define SRC_SCHEDULERS_HPP_

// C/C++
#include <unordered_map>
#include <utility>
#include <vector>

class MeshBlock;
class ParameterInput;

using IntegrationStage = std::pair<int, int>;
using TaskFunc = bool (*)(MeshBlock *, IntegrationStage);

struct TaskInfo {
  bool done = false;
  int load = 0;
  std::vector<TaskFunc> deps;
};

class Scheduler {
 public:
  /// constructor and destructor
  explicit Scheduler(MeshBlock *pmb);
  virtual ~Scheduler() {}

  /// functions
  bool DoTask(TaskFunc func);
  bool CheckDone(std::vector<TaskFunc> const &deps);
  void AddTasks(std::vector<TaskFunc> const &tasks);

 protected:
  std::unordered_map<TaskFunc, TaskInfo> tasks_;
  TaskFunc current_task_;
  IntegrationStage current_stage_;

 private:
  MeshBlock *pmy_block_;
};

using SchedulerPtr = std::shared_ptr<Scheduler>;

class SchedulerFactory {
 public:
  static SchedulerPtr Create(MeshBlock *pmb, ParameterInput *pin);
};

#endif  // SRC_SCHEDULERS_HPP_
