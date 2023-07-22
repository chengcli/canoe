#ifndef PHYSICS_HPP
#define PHYSICS_HPP

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"

// Forward declarations
class MeshBlock;
class ParameterInput;
enum class TaskStatus;
struct PhysicsPackage;

template <typename T>
class TaskManager;

//! \brief manages all physics package data and functions
class Physics {
 public:
  // data
  MeshBlock *pmy_block;
  TaskManager<Physics> *ptm;

  // functions
  Physics(MeshBlock *pmb, ParameterInput *pin);
  ~Physics();
  void Initialize(AthenaArray<Real> const &w);
  size_t RestartDataSizeInBytes();
  size_t DumpRestartData(char *pdst);
  size_t LoadRestartData(char *psrc);

  void ApplyPhysicsPackages(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                            Real time, Real dt);

  // package functions
  TaskStatus RelaxBotTemperature(AthenaArray<Real> &du,
                                 AthenaArray<Real> const &w, Real time,
                                 Real dt);
  TaskStatus RelaxBotVelocity(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                              Real time, Real dt);
  TaskStatus RelaxBotComposition(AthenaArray<Real> &du,
                                 AthenaArray<Real> const &w, Real time,
                                 Real dt);

  TaskStatus TopSpongeLayer(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                            Real time, Real dt);
  TaskStatus BotSpongeLayer(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                            Real time, Real dt);
  TaskStatus LeftSpongeLayer(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                             Real time, Real dt);
  TaskStatus RightSpongeLayer(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                              Real time, Real dt);

  TaskStatus TopCooling(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                        Real time, Real dt);
  TaskStatus BotHeating(AthenaArray<Real> &du, AthenaArray<Real> const &w,
                        Real time, Real dt);

 protected:
  // package list
  std::vector<PhysicsPackage> packages_;

  // bottom boundary condition
  AthenaArray<Real> hydro_bot_;
  AthenaArray<Real> tem_bot_, vel_bot_, com_bot_;

  // parameters for relax bottom temperature
  Real Tbot_, tau_Tbot_, tau_Ubot_, tau_Qbot_;

  // parameters for sponge layer
  Real tau_top_;   /**< top sponge layer time scale [s] */
  Real width_top_; /**< top sponge layer width [m] */

  Real tau_bot_;   /**< bottom sponge layer time scale [s] */
  Real width_bot_; /**< bottom sponge layer width [m] */

  Real tau_left_;   /**< left sponge layer time scale [s] */
  Real width_left_; /**< left sponge layer width [m] */

  Real tau_right_;   /**< right sponge layer time scale [s] */
  Real width_right_; /**< right sponge layer width [m] */

  // parameters for heating/cooling
  Real dTdt_;     /**< heating rate [K/s] */
  Real flux_top_; /**< top heating flux [W/m^2] */
  Real flux_bot_; /**< bot heating flux [W/m^2] */

  // parameters for friction
  Real friction_bot_tau_;
};

//! \brief task to do on a meshblock
struct PhysicsPackage {
  uint64_t id;        //!< encodes task using bit positions
  uint64_t dep;       //!< encodes dependencies to other tasks
  uint64_t conflict;  //!< encodes conflict tasks
  int cost;           //!< cost of this task
  bool load_balance;  //!< whether to do load balance

  // ptr to member function
  TaskStatus (Physics::*Function)(AthenaArray<Real> &,
                                  AthenaArray<Real> const &, Real, Real);
};

namespace PhysicsPackageNames {
const uint64_t NONE = 0LL;
const uint64_t FIX_BOT_TEMPERATURE = 1LL << 0;
const uint64_t FIX_BOT_VELOCITY = 1LL << 1;
const uint64_t FIX_BOT_COMPOSITION = 1LL << 2;
const uint64_t TOP_SPONGE_LAYER = 1LL << 3;
const uint64_t BOT_SPONGE_LAYER = 1LL << 4;
const uint64_t LFT_SPONGE_LAYER = 1LL << 5;
const uint64_t RHT_SPONGE_LAYER = 1LL << 6;
const uint64_t TOP_COOLING = 1LL << 7;
const uint64_t BOT_HEATING = 1LL << 8;
}  // namespace PhysicsPackageNames

#endif
