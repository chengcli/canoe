// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

// harp
#include "harp/radiation.hpp"
#include "harp/radiation_band.hpp"

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"
#include "snap/turbulence/turbulence_model.hpp"

// microphysics
#include "microphysics/microphysics.hpp"

// flask
#include "flask/chemistry.hpp"

// tracer
#include "tracer/tracer.hpp"

// astro
#include "astro/celestrial_body.hpp"

// exo3
#include "exo3/cubed_sphere.hpp"

// diagnostics
#include "diagnostics/diagnostics.hpp"

// forcing
#include "forcing/forcing.hpp"

// canoe
#include "impl.hpp"
#include "index_map.hpp"
#include "virtual_groups.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // decomposition
  pdec = std::make_shared<Decomposition>(pmb);

  // implicit methods
  phevi = std::make_shared<ImplicitSolver>(pmb, pin);

  // microphysics
  pmicro = std::make_shared<Microphysics>(pmb, pin);

  // radiation
  prad = std::make_shared<Radiation>(pmb, pin);

  // chemistry
  pchem = std::make_shared<Chemistry>(pmb, pin);

  // tracer
  ptracer = std::make_shared<Tracer>(pmb, pin);

  // turbulence
  pturb = TurbulenceFactory::Create(pmb, pin);

  // diagnostics
  all_diags = DiagnosticsFactory::CreateFrom(pmb, pin);

  // forcings
  all_forcings = ForcingFactory::CreateFrom(pmb, pin);

  // cubed sphere
  pexo3 = std::make_shared<CubedSphere>(pmb);

  // planet
  planet = PlanetFactory::CreateFrom(pmb, pin);
}

MeshBlock::Impl::~Impl() {}
