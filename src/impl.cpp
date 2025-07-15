// yaml
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/parameter_input.hpp>

// kintera
#include <kintera/kinetics/kinetics.hpp>

// harp
#include <harp/radiation/radiation.hpp>

// canoe
#include <configure.h>

// snap
#include <snap/eos/ideal_moist.hpp>
#include <snap/sedimentation/sedimentation.hpp>

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"

// astro
#include "astro/celestrial_body.hpp"

// exo3
#include "exo3/cubed_sphere.hpp"

// canoe
#include "impl.hpp"
#include "virtual_groups.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // decomposition
  pdec = std::make_shared<Decomposition>(pmb);

  // implicit methods
  phevi = std::make_shared<ImplicitSolver>(pmb, pin);

  // thermodynamcis
  auto yaml_file = pin->GetOrAddString("problem", "config_file", "");
  auto op_thermo = kintera::ThermoOptions::from_yaml(yaml_file);

  // microphysics
  auto op_kinet = kintera::KineticsOptions::from_yaml(yaml_file);
  pkinet = std::make_shared<kintera::KineticsImpl>(op_kinet);

  // radiation
  auto op_rad = harp::RadiationOptions::from_yaml(yaml_file);
  prad = std::make_shared<harp::RadiationImpl>(op_rad);

  // coordinate
  auto config = YAML::LoadFile(yaml_file);
  auto op_coord = snap::CoordinateOptions::from_yaml(config["geometry"]);

  // eos
  auto op_eos = snap::EquationOfStateOptions::from_yaml(
      config["dynamics"]["equation-of-state"]);
  op_eos.coord() = op_coord;
  op_eos.thermo() = op_thermo;
  peos = std::make_shared<snap::IdealMoistImpl>(op_eos);

  // sedimentation
  snap::SedHydroOptions op_sed;
  op_sed.sedvel() = snap::SedVelOptions::from_yaml(config);
  op_sed.eos() = op_eos;
  psed = std::make_shared<snap::SedHydroImpl>(op_sed);

  // cubed sphere
  pexo3 = std::make_shared<CubedSphere>(pmb);

  // planet
  planet = PlanetFactory::CreateFrom(pmb, pin);
}

MeshBlock::Impl::~Impl() {}
