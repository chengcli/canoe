//! \file turbulence_tasks.cpp
//! \brief declared in task_list/task_list.hpp

// C/C++ headers
#include <iostream>

// Athena++ headers
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <impl.hpp>

// tasklist
#include <tasklist/extra_tasks.hpp>

// snap
#include "turbulence_model.hpp"

TaskStatus ImplicitHydroTasks::CalculateTurbulenceFlux(MeshBlock *pmb,
                                                       int stage) {
  // std::cout << "calculate turbulence flux" << std::endl;
  auto pturb = pmb->pimpl->pturb;
  if (stage <= nstages) {
    pturb->calculateFluxes(pturb->r, 2);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::SendTurbulenceFlux(MeshBlock *pmb, int stage) {
  // std::cout << "send turbulence flux" << std::endl;
  pmb->pimpl->pturb->sbvar.SendFluxCorrection();
  return TaskStatus::success;
}

TaskStatus ImplicitHydroTasks::ReceiveTurbulenceFlux(MeshBlock *pmb,
                                                     int stage) {
  // std::cout << "receive turbulence flux" << std::endl;
  if (pmb->pimpl->pturb->sbvar.ReceiveFluxCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus ImplicitHydroTasks::IntegrateTurbulence(MeshBlock *pmb, int stage) {
  // std::cout << "integrate turbulence" << std::endl;
  auto pturb = pmb->pimpl->pturb;
  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to
    // IntegrateHydro, IntegrateField
    Real ave_wghts[5];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage - 1].delta;
    ave_wghts[2] = 0.0;
    ave_wghts[3] = 0.0;
    ave_wghts[4] = 0.0;
    pmb->WeightedAve(pturb->s1, pturb->s, pturb->s2, pturb->s2, pturb->s2,
                     ave_wghts);

    ave_wghts[0] = stage_wghts[stage - 1].gamma_1;
    ave_wghts[1] = stage_wghts[stage - 1].gamma_2;
    ave_wghts[2] = stage_wghts[stage - 1].gamma_3;
    if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0)
      pturb->s.SwapAthenaArray(pturb->s1);
    else
      pmb->WeightedAve(pturb->s, pturb->s1, pturb->s2, pturb->s2, pturb->s2,
                       ave_wghts);

    const Real wght = stage_wghts[stage - 1].beta * pmb->pmy_mesh->dt;
    pturb->addFluxDivergence(wght, pturb->s);

    // turbulence forcing
    Real dt = (stage_wghts[(stage - 1)].beta) * (pmb->pmy_mesh->dt);
    pturb->driveTurbulence(pturb->s, pturb->r, pmb->phydro->w, dt);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::SendTurbulence(MeshBlock *pmb, int stage) {
  // std::cout << "send turbulence" << std::endl;
  auto pturb = pmb->pimpl->pturb;
  if (stage <= nstages) {
    // Swap TurbulenceModel quantity in BoundaryVariable interface back to
    // conserved var formulation (also needed in SetBoundariesTurbulence() since
    // the tasks are independent)
    pturb->sbvar.var_cc = &(pturb->s);
    pturb->sbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

TaskStatus ImplicitHydroTasks::ReceiveTurbulence(MeshBlock *pmb, int stage) {
  // std::cout << "recv turbulence" << std::endl;
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pimpl->pturb->sbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

TaskStatus ImplicitHydroTasks::SetBoundariesTurbulence(MeshBlock *pmb,
                                                       int stage) {
  // std::cout << "set turbulence boundary" << std::endl;
  auto pturb = pmb->pimpl->pturb;
  if (stage <= nstages) {
    // Set Turbulence quantity in BoundaryVariable interface to cons var
    // formulation
    pturb->sbvar.var_cc = &(pturb->s);
    pturb->sbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}
