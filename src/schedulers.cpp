//! \file schedulers.cpp
//! \brief Implementation of scheduler

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include "schedulers.hpp"

Scheduler::Scheduler(MeshBlock* pmb) : pmy_block_(pmb) {}

bool Scheduler::DoTask(TaskFunc func) {
  current_task_ = func;
  return func(pmy_block_, current_stage_);
}

bool Scheduler::CheckDone(std::vector<TaskFunc> const& deps) {
  if (deps.empty()) return true;

  tasks_[current_task_].deps = deps;

  for (auto& dep_task : deps) {
    auto it = tasks_.find(dep_task);

    // find dependent task
    if (it != tasks_.end()) {
      if (it->second.done) continue;
    } else {  // register task, initialize to todo
      tasks_[dep_task] = {false, {}};
    }

    TaskFunc parent_task = current_task_;
    tasks_[dep_task].done = DoTask(dep_task);
    current_task_ = parent_task;

    if (!tasks_[dep_task].done) return false;
  }

  return true;
}

void Scheduler::AddTasks(std::vector<TaskFunc> const& tasks) {
  for (auto& task : tasks) {
    auto it = tasks_.find(task);

    // did not find task
    if (it == tasks_.end()) {
      tasks_[task] = {false, {}};
    }
  }
}

SchedulerPtr SchedulerFactory::Create(MeshBlock* pmb, ParameterInput* pin) {
  SchedulerPtr scheduler;
  scheduler = std::make_shared<Scheduler>(pmb);
  return scheduler;
}
