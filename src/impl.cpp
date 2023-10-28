// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

// application
#include <application/exceptions.hpp>

// harp
#include "harp/radiation.hpp"
#include "harp/radiation_band.hpp"

// inversion
#include "inversion/inversion.hpp"
#include "inversion/inversion_helper.hpp"

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"
#include "snap/turbulence/turbulence_model.hpp"

// microphysics
#include "microphysics/microphysics.hpp"

// n-body
#include "nbody/particles.hpp"

// canoe
#include "impl.hpp"
#include "index_map.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // decomposition
  pdec = std::make_shared<Decomposition>(pmb);

  // implicit methods
  phevi = std::make_shared<ImplicitSolver>(pmb, pin);

  // microphysics
  pmicro = std::make_shared<Microphysics>(pmb, pin);

  // chemistry
  pchem = std::make_shared<Chemistry>(pmb, pin);

  // tracer
  ptracer = std::make_shared<Tracer>(pmb, pin);

  // turbulence
  if (pin->DoesParameterExist("hydro", "turbulence")) {
    std::string turbulence_model = pin->GetString("hydro", "turbulence");
    if (turbulence_model == "none") {
      pturb = std::make_shared<TurbulenceModel>(pmb, pin);
    } else if (turbulence_model == "kepsilon") {
      pturb = std::make_shared<KEpsilonTurbulence>(pmb, pin);
      if (NTURBULENCE < 2)
        throw NotImplementedError(
            "NTURBULENCE must be at least 2 for k-epsilon model");
    } else {
      throw NotImplementedError(turbulence_model);
    }
  }

  // radiation
  prad = std::make_shared<Radiation>(pmb, pin);

  // inversion queue
  // all_fits = Inversion::NewInversionQueue(pmb, pin);

  // particles
  all_particles = ParticlesFactory::create_all_particles(pmb, pin);
}

MeshBlock::Impl::~Impl() {}

int find_pressure_level_lesser(Real pres, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie) {
  for (int i = is; i <= ie; ++i)
    if (w(IPR, k, j, i) < pres) return i;

  return ie + 1;
}
