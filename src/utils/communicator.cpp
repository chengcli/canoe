// C/C++
#include <iostream>
#include <mutex>

// athena
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// utils
#include "communicator.hpp"

static std::mutex comm_mutex;

Communicator *Communicator::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(comm_mutex);

  if (mycomm_ == nullptr) {
    mycomm_ = new Communicator();
  }

  return mycomm_;
}

Communicator *Communicator::mycomm_ = nullptr;
