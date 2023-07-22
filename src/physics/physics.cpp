// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../task_list/task_manager.hpp"
#include "physics.hpp"

Physics::~Physics() { delete ptm; }

//! \todo THIS HAS BEEN CHANGED u -> du. CHECK other packages for updates
void Physics::ApplyPhysicsPackages(AthenaArray<Real> &du,
                                   AthenaArray<Real> const &w, Real time,
                                   Real dt) {
  std::stringstream msg;
  int count = 0;
  ptm->Reset();

  while (count < 100) {
    TaskListStatus status = ptm->DoNextJob(du, w, time, dt, packages_);
    if (status == TaskListStatus::complete) break;
    count++;
  }

  if (count >= 100) {
    msg << "### FATAL ERROR in Physics::ApplyPhysicsPackages" << std::endl
        << "Physics Package stuck." << std::endl;
    ATHENA_ERROR(msg);
  }
}

size_t Physics::RestartDataSizeInBytes() {
  size_t size = 0;
  size += hydro_bot_.GetSizeInBytes();

  return size;
}

size_t Physics::DumpRestartData(char *pdst) {
  std::memcpy(pdst, hydro_bot_.data(), hydro_bot_.GetSizeInBytes());
  return RestartDataSizeInBytes();
}

size_t Physics::LoadRestartData(char *psrc) {
  std::memcpy(hydro_bot_.data(), psrc, hydro_bot_.GetSizeInBytes());
  return RestartDataSizeInBytes();
}
