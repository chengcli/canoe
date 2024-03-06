// C/C++
#include <array>
#include <cmath>

// external
#include <application/application.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/stride_iterator.hpp>

// canoe
#include <configure.hpp>

// climath
#include <climath/broyden_root.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// scm
#include "single_column.hpp"

struct TPBottomSolver {
  SingleColumn *pscm;
  int il, iu;
  Real mass0, enthalpy0;
  AirParcel air0;
};

//! wrapper function passing to broyden_root
//! over SingleColumn::findTPBottom
void find_tp_bottom(int n, double *x, double *f, void *arg) {
  TPBottomSolver *psolver = static_cast<TPBottomSolver *>(arg);

  AirParcel air = psolver->air0;

  air.w[IDN] *= x[0];
  air.w[IPR] *= x[1];

  auto result =
      psolver->pscm->findAdiabaticMassEnthalpy(air, psolver->il, psolver->iu);

  f[0] = result[0] / psolver->mass0 - 1.;
  f[1] = result[1] / psolver->enthalpy0 - 1.;
}

std::array<Real, 2> SingleColumn::findAdiabaticMassEnthalpy(AirParcel air,
                                                            int il, int iu) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmy_block_->pcoord;
  Real grav = -pmy_block_->phydro->hsrc.GetG1();

  // adjust pressure and density, and examine mass and energy conservation
  Real mass1 = 0., enthalpy1 = 0.;

  for (int i = il; i <= iu; ++i) {
    air.ToMassConcentration();

    Real density = air.w[IDN];
    for (int n = 1; n <= NVAPOR; ++n) {
      density += air.w[n];
    }

    mass1 += density * vol_(il);
    enthalpy1 += (air.w[IEN] + density * grav * pcoord->x1v(il)) * vol_(il);
    pthermo->Extrapolate(&air, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
  }

  return std::array<Real, 2>({mass1, enthalpy1});
}

std::array<Real, 2> SingleColumn::FindMassEnthalpy(AirColumn const &ac, int k,
                                                   int j, int il, int iu) {
  auto pcoord = pmy_block_->pcoord;
  Real grav = -pmy_block_->phydro->hsrc.GetG1();

  pcoord->CellVolume(k, j, il, iu, vol_);

  Real mass = 0.;
  Real enthalpy = 0.;

  for (int i = il; i <= iu; ++i) {
    auto air = ac[i];
    air.ToMassConcentration();

    Real density = air.w[IDN];
    for (int n = 1; n <= NVAPOR; ++n) {
      density += air.w[n];
    }
    mass += density * vol_(i);
    enthalpy += (air.w[IEN] + density * grav * pcoord->x1v(i)) * vol_(i);
  }

  return std::array<Real, 2>({mass, enthalpy});
}

// type of the function
void SingleColumn::FindUnstableRange(AirColumn const &ac, int il, int iu,
                                     std::vector<std::array<int, 2>> &ranges) {
  if (il >= iu) return;

  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmy_block_->pcoord;
  auto phydro = pmy_block_->phydro;

  Real den_tol = GetPar<Real>("den_tol");
  Real grav = -pmy_block_->phydro->hsrc.GetG1();

  // determine il where convective adjustment starts
  for (; il <= iu; ++il) {
    auto air0 = ac[il];
    auto air1 = ac[il + 1];
    if (air0.w[IDN] < 0. || air0.w[IPR] < 0.) break;
    pthermo->Extrapolate(&air0, pcoord->dx1f(il),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air0.ToMassFraction();
    air1.ToMassFraction();

    if (air0.w[IDN] - air1.w[IDN] < -den_tol) break;
  }

  // determine iu where convective adjustment stops
  auto air0 = ac[il];

  int i = il;
  for (; i < iu; ++i) {
    auto air1 = ac[i + 1];
    pthermo->Extrapolate(&air0, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
    air0.ToMassFraction();
    air1.ToMassFraction();

    if (air1.w[IDN] < 0. || air1.w[IPR] < 0.) continue;
    if (air0.w[IDN] - air1.w[IDN] > den_tol) break;
  }

  ranges.push_back(std::array<int, 2>({il, i}));

  FindUnstableRange(ac, i + 1, iu, ranges);
}

void SingleColumn::ConvectiveAdjustment(AirColumn &ac, int k, int j, int il,
                                        int iu) {
  if (il >= iu) return;

  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmy_block_->pcoord;
  Real grav = -pmy_block_->phydro->hsrc.GetG1();

  Real tp_tol = GetPar<Real>("rel_tol");
  int max_iter = GetPar<int>("max_iter");

  TPBottomSolver solver;
  solver.pscm = this;
  solver.il = il;
  solver.iu = iu;
  solver.air0 = ac[il];
  solver.air0.ToMassFraction();

  // sum the energy and mass of all air parcels that to be adjusted
  // (conservation of energy and mass)
  auto result = FindMassEnthalpy(ac, k, j, il, iu);
  solver.mass0 = result[0];
  solver.enthalpy0 = result[1];

  if (solver.mass0 < 0. || solver.enthalpy0 < 0.) {
    throw std::runtime_error("ConvectiveAdjustment: negative mass or enthalpy");
  }

  Real tp_bot[2] = {1., 1.};
  int status =
      broyden_root(2, tp_bot, find_tp_bottom, tp_tol, max_iter, &solver);
  /*if (status != 0) {
    if (status == 4) {
      std::cout << "status = " << status << std::endl;
      throw std::runtime_error(
          "ConvectiveAdjustment: "
          "Broyden's method failed to converge");
    } else {
      std::cout << "maximum iteration exceeded" << std::endl;
    }
  }*/

  solver.air0.w[IDN] *= tp_bot[0];
  solver.air0.w[IPR] *= tp_bot[1];

  auto app = Application::Logger("single_column");

  for (int i = il; i <= iu; ++i) {
    solver.air0.ConvertTo(ac[i].GetType());
    ac[i].w[IDN] = solver.air0.w[IDN];
    ac[i].w[IPR] = solver.air0.w[IPR];
    ac[i].w[IVX] /= 2.;
    ac[i].w[IVY] /= 2.;
    ac[i].w[IVZ] /= 2.;

    std::stringstream msg;
    msg << "ac[" << k << "," << j << "," << i << "] = " << ac[i];
    app->Log(msg.str());

    pthermo->Extrapolate(&solver.air0, pcoord->dx1f(i),
                         Thermodynamics::Method::ReversibleAdiabat, grav);
  }
}
