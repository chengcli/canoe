//! \file virtual_groups.cpp
//! \brief Implementation of virtual groups

// canoe
#include "virtual_groups.hpp"

std::unordered_map<TaskFunc, TaskInfo> TasksGroup::tasks_ = {};

bool TasksGroup::After(std::vector<TaskFunc> const& deps) {
  if (deps.empty()) return true;

  tasks_[current_task_].second = deps;

  for (auto& dep_task : deps) {
    auto it = tasks_.find(dep_task);

    // find dependent task
    if (it != tasks_.end()) {
      if (it->second.first) continue;
    } else {  // register task, initialize to todo
      tasks_[dep_task] = {false, {}};
    }

    TaskFunc parent_task = current_task_;
    current_task_ = dep_task;
    tasks_[dep_task].first = DoTask(dep_task);
    current_task_ = parent_task;

    if (!tasks_[dep_task].first) return false;
  }

  return true;
}

void TasksGroup::addTasks(std::vector<TaskFunc> const& tasks) {
  for (auto& task : tasks) {
    auto it = tasks_.find(task);

    // did not find task
    if (it == tasks_.end()) {
      tasks_[task] = {false, {}};
    }
  }
}
