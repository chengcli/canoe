// athena
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>

#include <random>

Real driving_acceleration;

template<typename T>
auto square(const T x) {
  return x * x;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1);
  SetUserOutputVariableName(0, "temp");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
      }
}

inline bool is_wall(MeshBlock *pmb, const int j) {
  const Real x2 = pmb->pcoord->x2v(j);
  const Real dx2 = pmb->pcoord->dx2f(j);
  return (
    (
      x2 < pmb->pmy_mesh->mesh_size.x2min + dx2
    ) || (
      x2 > pmb->pmy_mesh->mesh_size.x2max - dx2
    )
  );
}


void WallInteraction(MeshBlock *pmb, Real const time, Real const dt,
                     AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                     AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                     AthenaArray<Real> &s) {
  auto pthermo = Thermodynamics::GetInstance();

  const Real nu_iso = pmb->phydro->hdif.nu_iso;

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        if (is_wall(pmb, j)) {
          const auto w_kji = w.at(k, j, i);
          const Real dy = pmb->pcoord->dx2f(j);
          const Real r = -dt * nu_iso / square(0.5 * dy);

          const int nvs[] = {IVX, IVZ};
          for (auto n: nvs) {
            u(n, k, j, i) += r * w_kji[n] * w_kji[IDN];
          }
        }
      }
    }
  }
}

void DrivingAcceleration(MeshBlock *pmb, Real const time, Real const dt,
                     AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                     AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                     AthenaArray<Real> &s) {
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        u(IVX, k, j, i) += dt * driving_acceleration * w(IDN, k, j, i);
      }
    }
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  WallInteraction(pmb, time, dt, w, r, bcc, u, s);
  DrivingAcceleration(pmb, time, dt, w, r, bcc, u, s);
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  // index
  driving_acceleration = pin->GetReal("problem", "driving_acceleration");

  EnrollUserExplicitSourceFunction(Forcing);
}

template<class Real>
class IdealGas {
  public:
    Real gas_constant;
    Real specific_cv;

    IdealGas(const Real gas_constant, const Real specific_cv):
      gas_constant(gas_constant), specific_cv(specific_cv) {
    }

    template<class R1, class R2>
    inline auto density(const R1 &temp, const R2 &pres) const {
      return pres / (gas_constant * temp);
    }

    template<class R1>
    inline auto specific_internal_energy(const R1 &temp) const {
      return specific_cv * temp;
    }

    template<class R1>
    inline auto specific_enthalpy(const R1 &temp) const {
      return (specific_cv + gas_constant) * temp;
    }
};

template<class Real>
class CondensedMatter {
  public:
    IdealGas<Real> gas;
    Real temp3;
    Real pres3;
    Real beta;
    Real delta;

    CondensedMatter(const IdealGas<Real> &gas,
        const Real temp3, const Real pres3,
        const Real beta, const Real delta):
      gas(gas), temp3(temp3), pres3(pres3), beta(beta), delta(delta) {}

    template<class R1>
    inline auto specific_internal_energy(const R1 &temp) const {
      return (
        gas.specific_enthalpy(temp)
        + gas.gas_constant * (-beta * temp3 + delta * temp)
      );
    }

    template<class R>
    inline auto pres_sat(const R &temp) const {
      auto t3 = temp / temp3;
      return pres3 * exp(beta * (1. - 1./t3) - delta * log(t3));
    }

    template<class R>
    inline auto vapor_density_sat(const R &temp) const {
      return gas.density(temp, pres_sat(temp));
    }
};

inline auto WaterIceEOS() {
  const double Avogadro = 6.02214076e23;
  const double Boltzmann = 1.380649e-23;
  const Real atomic_mass_H = 1.008e-3;
  const Real atomic_mass_O = 15.999e-3;
  const Real universial_gas_constant = Avogadro * Boltzmann;

  const Real water_mw = 2 * atomic_mass_H  + atomic_mass_O;

  const Real water_gas_constant = universial_gas_constant / water_mw;

  const Real water_vapor_cp_mol = 37.4;

  const Real water_vapor_cp = water_vapor_cp_mol / water_mw;

  const Real water_vapor_cv = water_vapor_cp - water_gas_constant;

  const Real temp3 = 273.16;
  const Real pres3 = 611.7;
  const Real beta = 24.845;
  const Real delta = 4.986;

  IdealGas<Real> water_vapor(water_gas_constant, water_vapor_cv);
  CondensedMatter<Real> water_ice(water_vapor, temp3, pres3, beta, delta);
  return water_ice;
}


void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  const auto mesh_size = pmy_mesh->mesh_size;
  const Real yc = 0.5 * (mesh_size.x2max + mesh_size.x2min);
  const Real yd = 0.5 * (mesh_size.x2max - mesh_size.x2min);

  auto pthermo = Thermodynamics::GetInstance();
  auto water_ice_eos = WaterIceEOS();
  const Real ice_fraction = pin->GetReal("initialcondition", "ice_fraction");
  const Real temperature = pin->GetReal("initialcondition", "temperature");
  const Real pressure = water_ice_eos.pres_sat(temperature);
  const Real density = (
    water_ice_eos.vapor_density_sat(temperature)
    / (1. - ice_fraction)
  );

  const Real nu_iso = pin->GetReal("problem", "nu_iso");

  const Real uc = 0.5 * square(yd) * driving_acceleration / nu_iso;

  std::mt19937 mt(1234);
  std::uniform_real_distribution<Real> phi(0., 2 * M_PI);

  // populate to 3D mesh
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN, k, j, i) = density;
        phydro->w(pthermo->SpeciesIndex("H2O"), k, j, i) = 1. - ice_fraction;
        phydro->w(pthermo->SpeciesIndex("H2O(s)"), k, j, i) = ice_fraction;
        phydro->w(IVX, k, j, i) = uc * (1. - square((pcoord->x2v(j) - yc) / yd));
        phydro->w(IPR, k, j, i) = pressure;
      }
    }
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
