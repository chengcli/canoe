#include <iostream>
#include <vector>

#include "ice_shell.hpp"

int main(int argc, char* argv[]) {

  typedef double Real;

  const int nx = 10;
  const int nz = 10;
  const int nb = nx + nz;

  const Real vapor_p3 = 611.7; // Pascal
  const Real vapor_t3 = 273.16; // Kelvin
  const Real vapor_gas_const = 461.5; // J kg-1 K-1
  const Real vapor_adiabatic_index = 1.4;
  const Real vapor_beta = 24.845;
  const Real vapor_delta = 4.986;

  const Real outer_surface_temp_eff = 67.; // Kelvin
  const Real stefan_boltzmann_const = 5.67e-8;

  const Real ocean_temp = 273.16; // Kelvin
  const Real ice_k = 651.; // Watts per meter

  IceShell::VaporCondensation<Real> cond
    (vapor_p3, vapor_t3, vapor_gas_const, vapor_adiabatic_index,
       vapor_beta, vapor_delta);
  IceShell::OuterSurfaceRadiation<Real> rad
    (outer_surface_temp_eff, stefan_boltzmann_const);
  IceShell::TemperatureDependentDiffusivity<Real> tdd
    (ocean_temp, ice_k);

  // random exponential grid
  std::vector<Real> dx;
  std::vector<Real> dz;
  Real d = 1.;

  for(int i = 0; i < nx; ++i) {
    dx.push_back(d);
    d *= 1.05;
  }

  for(int k = 0; k < nz; ++k) {
    dz.push_back(d);
    d *= 1.05;
  }

  IceShell::IceBoundaryModel<Real> ice_boundary_model
    (dx, dz, cond, rad, tdd);

  IceShell::BoundaryValue<Real> air_temp (nx, nz);
  IceShell::BoundaryValue<Real> ice_temp (nx, nz);


  // random air temperature
  for (int i = 0; i < nb; ++i) {
    Real t = std::max(270. - i * 200./(nx + nz), 160.);
    if (i < nz) {
      air_temp.set_side(i, t);
    } else {
      air_temp.set_top(i - nz, t);
    }
  }

  ice_boundary_model.solve(ice_temp, air_temp);

  std::cout << air_temp.data << std::endl;
  std::cout << ice_temp.data << std::endl;

  return 0;
}
