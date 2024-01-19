#ifndef SRC_TASKLIST_TASK_LIST_FACTORY_HPP_
#define SRC_TASKLIST_TASK_LIST_FACTORY_HPP_

// C/C++
#include <memory>

// application
#include <application/exceptions.hpp>

#include "extra_tasks.hpp"

using TaskListPtr = std::unique_ptr<TaskList>;

class TaskListFactory {
 public:
  static TaskListPtr CreateFrom(ParameterInput* pin, Mesh* mesh) {
    std::string tasklist_name;
    TaskListPtr ptlist;

    if (pin->DoesParameterExist("job", "tasklist")) {
      tasklist_name = pin->GetString("job", "tasklist");
    } else {
      tasklist_name = "ImplicitHydroTasks";
    }

    try {
      if (tasklist_name == "TimeIntegratorTaskList") {
        ptlist = std::make_unique<TimeIntegratorTaskList>(pin, mesh);
      } else if (tasklist_name == "ImplicitHydroTasks") {
        ptlist = std::make_unique<ImplicitHydroTasks>(pin, mesh);
      } else if (tasklist_name == "InversionTasks") {
        ptlist = std::make_unique<InversionTasks>(pin, mesh);
      } else {
        throw RuntimeError("main", "Unknown tasklist name: " + tasklist_name);
      }
    } catch (std::bad_alloc& ba) {
      throw RuntimeError("main", "TaskList memory allocation failed");
    }

    return ptlist;
  }
};

#endif  // SRC_TASKLIST_TASK_LIST_FACTORY_HPP_
