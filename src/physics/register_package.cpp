// Athena++ headers
#include "../athena.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../task_list/task_manager.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "physics.hpp"

using namespace PhysicsPackageNames;

Physics::Physics(MeshBlock *pmb, ParameterInput *pin) : pmy_block(pmb) {
  pmb->pdebug->Enter("Physics");
  std::stringstream &msg = pmb->pdebug->msg;
  char package_names[1024], *p;
  std::string str = pin->GetOrAddString("physics", "packages", "");
  std::strcpy(package_names, str.c_str());
  p = std::strtok(package_names, " ,");

  ptm = new TaskManager<Physics>(this);

  PhysicsPackage pkg;
  while (p != NULL) {
    if (std::strcmp(p, "fix_bot_temperature") == 0) {
      msg << "- use physcis fix_bot_temperature" << std::endl;
      pkg.id = FIX_BOT_TEMPERATURE;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::RelaxBotTemperature;
      ptm->AddPackage(pkg, "fix_bot_temperature");

      tau_Tbot_ = pin->GetReal("physics", "fix_bot_temperature.tau");
      // -1 means to use the initial condition
      Tbot_ = pin->GetOrAddReal("physics", "bot_temperature", -1);
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      msg << "- tau = " << tau_Tbot_ << " s" << std::endl;
      tem_bot_.InitWithShallowSlice(hydro_bot_, 3, IDN, 1);
    } else if (std::strcmp(p, "fix_bot_velocity") == 0) {
      msg << "- use physcis fix_bot_velocity" << std::endl;
      pkg.id = FIX_BOT_VELOCITY;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      // pkg.Function = &Physics::RelaxBotVelocity;
      ptm->AddPackage(pkg, "fix_bot_velocity");

      tau_Ubot_ = pin->GetReal("physics", "fix_bot_velocity.tau");
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      vel_bot_.InitWithShallowSlice(hydro_bot_, 3, IVX, 3);
    } else if (std::strcmp(p, "fix_bot_composition") == 0) {
      msg << "- use physcis fix_bot_composition" << std::endl;
      pkg.id = FIX_BOT_COMPOSITION;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      // pkg.Function = &Physics::RelaxBotComposition;
      ptm->AddPackage(pkg, "fix_bot_composition");

      tau_Ubot_ = pin->GetReal("physics", "fix_bot_composition.tau");
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      com_bot_.InitWithShallowSlice(hydro_bot_, 3, IDN, 1 + NVAPOR);
    } else if (std::strcmp(p, "top_sponge") == 0) {
      msg << "- use physcis top_sponge" << std::endl;
      pkg.id = TOP_SPONGE_LAYER;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::TopSpongeLayer;
      ptm->AddPackage(pkg, "top_sponge");

      tau_top_ = pin->GetReal("physics", "top_sponge.tau");
      width_top_ = pin->GetReal("physics", "top_sponge.width");
    } else if (std::strcmp(p, "bot_sponge") == 0) {
      msg << "- use physcis bot_sponge" << std::endl;
      pkg.id = BOT_SPONGE_LAYER;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::BotSpongeLayer;
      ptm->AddPackage(pkg, "bot_sponge");

      tau_bot_ = pin->GetReal("physics", "bot_sponge.tau");
      width_bot_ = pin->GetReal("physics", "bot_sponge.width");
    } else if (std::strcmp(p, "left_sponge") == 0) {
      msg << "- use physcis left_sponge" << std::endl;
      pkg.id = LFT_SPONGE_LAYER;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::LeftSpongeLayer;
      ptm->AddPackage(pkg, "left_sponge");

      tau_left_ = pin->GetReal("physics", "left_sponge.tau");
      width_left_ = pin->GetReal("physics", "left_sponge.width");
    } else if (std::strcmp(p, "right_sponge") == 0) {
      msg << "- use physcis right_sponge" << std::endl;
      pkg.id = RHT_SPONGE_LAYER;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::RightSpongeLayer;
      ptm->AddPackage(pkg, "right_sponge");

      tau_right_ = pin->GetReal("physics", "right_sponge.tau");
      width_right_ = pin->GetReal("physics", "right_sponge.width");
    } else if (std::strcmp(p, "top_cooling") == 0) {
      msg << "- use physcis top_cooling" << std::endl;
      pkg.id = TOP_COOLING;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::TopCooling;
      ptm->AddPackage(pkg, "top_cooling");

      // dTdt_ = pin->GetReal("physics", "top_cooling.rate")/86400.; // K/day to
      // K/s
      flux_top_ = pin->GetReal("physics", "top_cooling.flux");
    } else if (std::strcmp(p, "bot_heating") == 0) {
      msg << "- use physcis bot_heating" << std::endl;
      pkg.id = BOT_HEATING;
      pkg.dep = 0LL;
      pkg.conflict = 0LL;
      pkg.Function = &Physics::BotHeating;
      ptm->AddPackage(pkg, "bot_heating");

      flux_bot_ = pin->GetReal("physics", "bot_heating.flux");
    } else {
      msg << "### FATAL ERROR in function Physics::Physics" << std::endl
          << "Package '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    packages_.push_back(pkg);
    p = std::strtok(NULL, " ,");
  }
  pmb->pdebug->Leave();
}

void Physics::Initialize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block;

  /* find top and bot neighbor
  NeighborBlock ntop, nbot;
  bool has_top_neighbor = false;
  bool has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      nbot = nb;
      has_bot_neighbor = true;
    } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      ntop = nb;
      has_top_neighbor = true;
    }
  }

  if (has_bot_neighbor)
    ptm->RemoveTask(FIX_BOT_TEMPERATURE | FIX_BOT_VELOCITY | FIX_BOT_COMPOSITION
                  | BOT_HEATING);

  if (has_top_neighbor)
    ptm->RemoveTask(TOP_COOLING);*/

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      if (ptm->HasTask(FIX_BOT_TEMPERATURE))
        tem_bot_(k, j) = pmb->pthermo->GetTemp(w.at(k, j, pmb->is));
      if (ptm->HasTask(FIX_BOT_VELOCITY)) {
        vel_bot_(0, k, j) = w(IVX, k, j, pmb->is);
        vel_bot_(1, k, j) = w(IVY, k, j, pmb->is);
        vel_bot_(2, k, j) = w(IVZ, k, j, pmb->is);
      }
      if (ptm->HasTask(FIX_BOT_COMPOSITION))
        for (int n = 1; n <= NVAPOR; ++n)
          com_bot_(n, k, j) = w(n, k, j, pmb->is);
    }
}
