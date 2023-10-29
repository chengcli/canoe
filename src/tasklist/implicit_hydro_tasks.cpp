// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>

// snap
#include <snap/implicit/implicit_solver.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// microphysics
#include <microphysics/microphysics.hpp>

// tasklist
#include "extra_tasks.hpp"

using TaskFunction = TaskStatus (TaskList::*)(MeshBlock *, int);

int find_task(Task const *task_list, int ntasks, TaskID const &id) {
  for (int i = 0; i < ntasks; ++i) {
    if (task_list[i].task_id == id) return i;
  }
  throw NotFoundError("find_task", "Task Function");
}

ImplicitHydroTasks::ImplicitHydroTasks(ParameterInput *pin, Mesh *pm)
    : TimeIntegratorTaskList(pin, pm) {
  using namespace HydroIntegratorTaskNames;  // NOLINT (build/namespace)
  int itask;

  // update IntegrateHydro
  itask = find_task(task_list_, ntasks, INT_HYD);
  task_list_[itask].TaskFunc =
      static_cast<TaskFunction>(&ImplicitHydroTasks::IntegrateHydro);

  // update AddSourceTerms
  itask = find_task(task_list_, ntasks, SRC_TERM);
  task_list_[itask].TaskFunc =
      static_cast<TaskFunction>(&ImplicitHydroTasks::AddSourceTerms);

  // inject canoe tasks
  // **ATHENA TASKS**: ... -> (CALC_HYDFLX | CALC_SCLRFLX) ->
  // **CANOE TASKS** : ADD_FLX_CONS ->
  // **ATHENA TASKS**: (INT_HYD | INT_SCLR) -> SRC_TERM ->
  // **CANOE TASKS** : IMPLICIT_CORR -> UPDATE_ALLCONS ->
  // **ATHENA TASKS**: (SEND_HYD | SEND_SCLR) -> ...
  if (NSCALARS > 0) {
    AddTask(ADD_FLX_CONS, (CALC_HYDFLX | CALC_SCLRFLX));
  } else {
    AddTask(ADD_FLX_CONS, CALC_HYDFLX);
  }

  AddTask(IMPLICIT_CORR, SRC_TERM);
  AddTask(UPDATE_ALLCONS, IMPLICIT_CORR);

  // update dependency
  if (ORBITAL_ADVECTION) {
    itask = find_task(task_list_, ntasks, SEND_HYDORB);
    task_list_[itask].dependency = UPDATE_ALLCONS;
  } else {
    itask = find_task(task_list_, ntasks, SEND_HYD);
    task_list_[itask].dependency = UPDATE_ALLCONS;

    itask = find_task(task_list_, ntasks, SETB_HYD);
    task_list_[itask].dependency = (RECV_HYD | UPDATE_ALLCONS);
  }

  if (pm->multilevel || SHEAR_PERIODIC) {  // SMR or AMR or shear periodic
    itask = find_task(task_list_, ntasks, SEND_HYDFLX);
    task_list_[itask].dependency = ADD_FLX_CONS;

    itask = find_task(task_list_, ntasks, RECV_HYDFLX);
    task_list_[itask].dependency = ADD_FLX_CONS;
  } else {
    itask = find_task(task_list_, ntasks, INT_HYD);
    task_list_[itask].dependency = ADD_FLX_CONS;
  }

  if (NSCALARS > 0) {
    if (!ORBITAL_ADVECTION) {
      itask = find_task(task_list_, ntasks, SEND_SCLR);
      task_list_[itask].dependency = UPDATE_ALLCONS;

      itask = find_task(task_list_, ntasks, SETB_SCLR);
      task_list_[itask].dependency = (RECV_SCLR | UPDATE_ALLCONS);
    }
  }
}

void ImplicitHydroTasks::AddTask(TaskID const &id, TaskID const &dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  using namespace HydroIntegratorTaskNames;

  if (id == ADD_FLX_CONS) {
    task_list_[ntasks].TaskFunc =
        static_cast<TaskFunction>(&ImplicitHydroTasks::AddFluxToConserved);
    task_list_[ntasks].lb_time = true;
  } else if (id == IMPLICIT_CORR) {
    task_list_[ntasks].TaskFunc =
        static_cast<TaskFunction>(&ImplicitHydroTasks::ImplicitCorrection);
    task_list_[ntasks].lb_time = true;
  } else if (id == UPDATE_ALLCONS) {
    task_list_[ntasks].TaskFunc =
        static_cast<TaskFunction>(&ImplicitHydroTasks::UpdateAllConserved);
    task_list_[ntasks].lb_time = true;
  } else {
    throw NotFoundError("AddTask", "Task ID");
  }

  ntasks++;
}

TaskStatus ImplicitHydroTasks::IntegrateHydro(MeshBlock *pmb, int stage) {
  auto ph = pmb->phydro;
  auto pf = pmb->pfield;
  auto phevi = pmb->pimpl->phevi;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve)
    return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage - 1].main_stage) {
      // This time-integrator-specific averaging operation logic is identical to
      // FieldInt
      Real ave_wghts[5];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage - 1].delta;
      ave_wghts[2] = 0.0;
      ave_wghts[3] = 0.0;
      ave_wghts[4] = 0.0;
      pmb->WeightedAve(ph->u1, ph->u, ph->u2, ph->u0, ph->fl_div, ave_wghts);

      ave_wghts[0] = stage_wghts[stage - 1].gamma_1;
      ave_wghts[1] = stage_wghts[stage - 1].gamma_2;
      ave_wghts[2] = stage_wghts[stage - 1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0)
        ph->u.SwapAthenaArray(ph->u1);
      else
        pmb->WeightedAve(ph->u, ph->u1, ph->u2, ph->u0, ph->fl_div, ave_wghts);

      const Real wght = stage_wghts[stage - 1].beta * pmb->pmy_mesh->dt;
      pmb->pimpl->du.ZeroClear();
      ph->AddFluxDivergence(wght, pmb->pimpl->du);
      // add coordinate (geometric) source terms
      pmb->pcoord->AddCoordTermsDivergence(wght, ph->flux, ph->w, pf->bcc,
                                           pmb->pimpl->du);

      // Hardcode an additional flux divergence weighted average for the
      // penultimate stage of SSPRK(5,4) since it cannot be expressed in a 3S*
      // framework
      if (stage == 4 && integrator == "ssprk5_4") {
        // From Gottlieb (2009), u^(n+1) partial calculation
        ave_wghts[0] = -1.0;  // -u^(n) coeff.
        ave_wghts[1] = 0.0;
        ave_wghts[2] = 0.0;
        const Real beta = 0.063692468666290;  // F(u^(3)) coeff.
        const Real wght_ssp = beta * pmb->pmy_mesh->dt;
        // writing out to u2 register
        pmb->WeightedAve(ph->u2, ph->u1, ph->u2, ph->u0, ph->fl_div, ave_wghts);
        ph->AddFluxDivergence(wght_ssp, ph->u2);
        // add coordinate (geometric) source terms
        pmb->pcoord->AddCoordTermsDivergence(wght_ssp, ph->flux, ph->w, pf->bcc,
                                             ph->u2);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::AddSourceTerms(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  PassiveScalars *ps = pmb->pscalars;

  // return if there are no source terms to be added
  if (!(ph->hsrc.hydro_sourceterms_defined) ||
      pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve)
    return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage - 1].main_stage) {
      // Time at beginning of stage for u()
      Real t_start_stage = pmb->pmy_mesh->time +
                           stage_wghts[(stage - 1)].sbeta * pmb->pmy_mesh->dt;
      // Scaled coefficient for RHS update
      Real dt = (stage_wghts[(stage - 1)].beta) * (pmb->pmy_mesh->dt);
      // Evaluate the source terms at the time at the beginning of the stage
      ph->hsrc.AddSourceTerms(t_start_stage, dt, ph->flux, ph->w, ps->r,
                              pf->bcc, pmb->pimpl->du, ps->s);
      // pmb->pphy->ApplyPhysicsPackages(phevi->du, ph->w, pmb->pmy_mesh->time,
      // pmb->pmy_mesh->dt);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::AddFluxToConserved(MeshBlock *pmb, int stage) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto prad = pmb->pimpl->prad;

  if (stage <= nstages) {
    if (stage_wghts[stage - 1].main_stage) {
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) {
          prad->AddRadiativeFlux(pmb->phydro, k, j, is, ie + 1);
        }

      if (stage == nstages) {           // last stage
        if (prad->GetCounter() > 0.) {  // reuse previous flux
          prad->DecrementCounter(pmb->pmy_mesh->dt);
        } else {  // update radiative flux
          prad->ResetCounter();
          for (int k = ks; k <= ke; ++k)
            for (int j = js; j <= je; ++j)
              prad->CalRadiativeFlux(pmb->pmy_mesh->time, k, j, is, ie + 1);
        }
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::ImplicitCorrection(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  auto phevi = pmb->pimpl->phevi;
  Real dt = pmb->pmy_mesh->dt;

  if (stage <= nstages) {
    if (stage_wghts[stage - 1].main_stage) {
      // do implicit coorection at every stage
      // X3DIR
      if ((phevi->implicit_flag & (1 << 2)) && (pmb->ncells3 > 1)) {
        phevi->SetDirection(X3DIR);
        if (phevi->implicit_flag & (1 << 3))
          phevi->FullCorrection(pmb->pimpl->du, ph->w,
                                stage_wghts[stage - 1].beta * dt);
        else
          phevi->PartialCorrection(pmb->pimpl->du, ph->w,
                                   stage_wghts[stage - 1].beta * dt);
      }

      // X2DIR
      if ((phevi->implicit_flag & (1 << 1)) && (pmb->ncells2 > 1)) {
        phevi->SetDirection(X2DIR);
        if (phevi->implicit_flag & (1 << 3))
          phevi->FullCorrection(pmb->pimpl->du, ph->w,
                                stage_wghts[stage - 1].beta * dt);
        else
          phevi->PartialCorrection(pmb->pimpl->du, ph->w,
                                   stage_wghts[stage - 1].beta * dt);
      }

      // X1DIR
      if (phevi->implicit_flag & 1) {
        phevi->SetDirection(X1DIR);
        if (phevi->implicit_flag & (1 << 3))
          phevi->FullCorrection(pmb->pimpl->du, ph->w,
                                stage_wghts[stage - 1].beta * dt);
        else
          phevi->PartialCorrection(pmb->pimpl->du, ph->w,
                                   stage_wghts[stage - 1].beta * dt);
      }

      Real wghts[5];
      wghts[0] = 1.;
      wghts[1] = 1.;
      wghts[2] = 0.;
      wghts[3] = 0.;
      wghts[4] = 0.;
      pmb->WeightedAve(ph->u, pmb->pimpl->du, ph->u2, ph->u2, ph->u2, wghts);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus ImplicitHydroTasks::UpdateAllConserved(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pimpl->pmicro->SetSedimentationVelocityFromConserved(pmb->phydro);
  } else {
    return TaskStatus::fail;
  }

  // only do chemistry and thermodynamcis at last rk step
  if (stage != nstages) return TaskStatus::next;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto pthermo = Thermodynamics::GetInstance();

  auto pmicro = pmb->pimpl->pmicro;

  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++) {
      auto &&ac = AirParcelHelper::gather_from_conserved(pmb, k, j, is, ie);

      // pmicro->AddFrictionalHeating(air_column);

      pmicro->EvolveSystems(ac, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

      pthermo->SaturationAdjustment(ac);

      AirParcelHelper::distribute_to_conserved(pmb, k, j, is, ie, ac);
    }

  return TaskStatus::success;
}
