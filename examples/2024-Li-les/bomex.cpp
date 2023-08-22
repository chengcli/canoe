/* -------------------------------------------------------------------------------------
 * SNAP-Earth Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 * Zhaowyi Shen, California Institute of Technology
 * Zhihong Tan, Princeton University, GFDL
 *
 * Year: 2022
 * Contact: chengcli@umich.edu
 * Reference: Siebesma et al., 2003 (BOMEX)
 * -------------------------------------------------------------------------------------
 */

// @sect3{Include files}

// These input files are just like those in the @ref 2d.straka, so additional
// comments are not required.
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../math/core.h"  // sqr
#include "../math/root.h"  // sqr
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../utils/utils.hpp"

// @sect3{Preamble}

// We need two global variables here, reference pressure and gravity
Real p0, grav;

// water vapor id
enum { iH2O = 1 };

// Same as that in @ref 2d.straka, make outputs of temperature and potential
// temperature.
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp", "temperature", "K");
  SetUserOutputVariableName(1, "theta", "potential temperature", "K");
}

// Set temperature and potential temperature.
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real gamma = peos->GetGamma();
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(phydro->w.at(k, j, i));
        user_out_var(1, k, j, i) =
            PotentialTemp(phydro->w.at(k, j, i), p0, pthermo);
      }
}

// BOMEX Large scale forcing
void BomexForcing(MeshBlock *pmb, Real const time, Real const dt,
                  AthenaArray<Real> const &w, AthenaArray<Real> const &bcc,
                  AthenaArray<Real> &u) {
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real vstar = 0.28;
      Real cv = pthermo->getSpecificCv(w.at(k, j, is));

      // surface sensible heat flux, Eq B1
      u(IEN, k, j, is) += dt * 8.E-3 * w(IDN, k, j, is) * cv / pcoord->dx1f(is);

      // surface latent heat flux, Eq B1
      u(IEN, k, j, is) += dt * 5.2E-5 * w(IDN, k, j, is) / pcoord->dx1f(is);

      // surface momentum flux, Eq B2
      Real um =
          sqr(w(IVX, k, j, is)) + sqr(w(IVY, k, j, is)) + sqr(w(IVZ, k, j, is));
      for (int n = IVX; n <= IVZ; ++n) {
        u(n, k, j, is) +=
            -dt * sqr(vstar) / um * w(n, k, j, is) / pcoord->dx1f(is);
      }

      for (int i = is; i <= ie; ++i) {
        cv = pthermo->getSpecificCv(w.at(k, j, i));

        // large scale cooling
        Real Qr;  // K/s
        if (pcoord->x1v(i) < 1500.) {
          Qr = -2.0 / 86400.;
        } else if (pcoord->x1v(i) < 3000.) {
          Qr = -2.0 / 86400. * (3000. - pcoord->x1v(i)) / 1500.;
        }
        u(IEN, k, j, i) -= w(IDN, k, j, i) * cv * Qr * dt;

        // large scale drying
        Real dqtdt;  // 1/s
        if (pcoord->x1v(i) < 300.) {
          dqtdt = -1.2E-8;
        } else if (pcoord->x1v(i) < 500.) {
          dqtdt = -1.2E-8 * (500. - pcoord->x1v(i)) / 200.;
        }
        u(iH2O, k, j, i) -= w(IDN, k, j, i) * dqtdt * dt;

        // large scale subsidence
        Real wsub;  // m/s
        if (pcoord->x1v(i) < 1500.) {
          wsub = -0.65E-2 * pcoord->x1v(i) / 1500.;
        } else if (pcoord->x1v(i) < 2100.) {
          wsub = -0.65E-2 * (2100. - pcoord->x1v(i)) / 600.;
        }

        // 1. vapor
        u(iH2O, k, j, i) -= w(IDN, k, j, i) * wsub *
                            (w(iH2O, k, j, i + 1) - w(iH2O, k, j, i - 1)) /
                            (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));

        // 2. momentum
        for (int n = IVX; n <= IVZ; ++n)
          u(n, k, j, i) -= w(IDN, k, j, i) * wsub *
                           (w(n, k, j, i + 1) - w(n, k, j, i - 1)) /
                           (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));

        // 3. energy
        Real mse1 = MoistStaticEnergy(w.at(k, j, i + 1),
                                      grav * pcoord->x1v(i + 1), pthermo);
        Real ke1 = 0.5 * (sqr(w(IVX, k, j, i + 1)) + sqr(w(IVY, k, j, i + 1)) +
                          sqr(w(IVZ, k, j, i + 1)));
        Real mse2 = MoistStaticEnergy(w.at(k, j, i - 1),
                                      grav * pcoord->x1v(i - 1), pthermo);
        Real ke2 = 0.5 * (sqr(w(IVX, k, j, i - 1)) + sqr(w(IVY, k, j, i - 1)) +
                          sqr(w(IVZ, k, j, i - 1)));
        u(IEN, k, j, i) -= w(IDN, k, j, i) * wsub * (mse1 + ke1 - mse2 - ke2) /
                           (pcoord->x1v(i + 1) - pcoord->x1v(i - 1));

        // geostrophic wind forcing
        Real ug = -10. + 1.8E-3 * pcoord->x1v(i);
        Real vg = 0.;
        Real fcor = 0.376E-4;
        u(IM2, k, j, i) += fcor * (w(IM3, k, j, i) - vg);
        u(IM3, k, j, i) += -fcor * (w(IM2, k, j, i) - ug);
      }
    }
}

// Initialize surface pressure from input file.
void Mesh::InitUserMeshData(ParameterInput *pin) {
  p0 = pin->GetReal("problem", "p0");
  grav = -pin->GetReal("hydro", "grav_acc1");
  EnrollUserExplicitSourceFunction(BomexForcing);
}

// theta_l Solver
struct ThetaLSolver {
  Thermodynamics *pthermo;
  Real **w2;
  Real dz;
  Real thetal;
  Real lv_ov_cpd;
  Real qliq;
};

Real solveDensityGivenThetaL(Real rho, void *aux) {
  // grav parameter is not used in hydrostatic formulation, set to zero
  ThetaLSolver *psolver = static_cast<ThetaLSolver *>(aux);
  Real **w2 = psolver->w2;
  Real lv_ov_cpd = psolver->lv_ov_cpd;
  Real qliq = psolver->qliq;
  Thermodynamics *pthermo = psolver->pthermo;
  w2[1][IPR] = w2[0][IPR] - 0.5 * (w2[0][IDN] + rho) * grav * psolver->dz;
  w2[1][IDN] = rho;
  //! no liquid water initially
  Real thetal = LiquidWaterPotentialTemp(w2[1], p0, lv_ov_cpd, qliq, pthermo);
  return thetal - psolver->thetal;
}

// @sect3{Initial condition}

// We do not need forcings other than gravity in this problem,
// so we go directly to the initial condition.
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // Similar to @ref 2d.straka, read variables in the input file
  Real gamma = pin->GetReal("hydro", "gamma");
  Real grav = -phydro->hsrc.GetG1();
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;
  Real Ts = pin->GetOrAddReal("problem", "Ts", 300.4);
  Real Ps = pin->GetOrAddReal("problem", "Ps", 1.015E5);
  Real qs = pin->GetOrAddReal("problem", "qs", 17.0) / 1.E3;

  Real **w2, qt;
  NewCArray(w2, 2, NHYDRO + 2 * NVAPOR);
  memset(*w2, 0., 2 * (NHYDRO + 2 * NVAPOR) * sizeof(Real));

  ThetaLSolver solver;
  solver.pthermo = pthermo;
  solver.w2 = w2;
  solver.lv_ov_cpd = 0.;
  solver.qliq = 0.;

  // Loop over the grids and set initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      // bottom boundary
      w2[0][iH2O] = qs;
      w2[0][IPR] = Ps;
      w2[0][IDN] = Ps / (Rd * pthermo->RovRd(w2[0]) * Ts);
      Real rho = w2[0][IDN];

      for (int i = is; i <= ie; ++i) {
        // solve for layer interface
        Real x1 = pcoord->x1f(i + 1);
        solver.dz = pcoord->dx1f(i);

        if (x1 < 520.) {
          qt = 17.0 * (520. - x1) / 520. + 16.3 * x1 / 520.;
          solver.thetal = 298.7;
        } else if (x1 < 1480.) {
          qt = 16.3 * (1480. - x1) / (1480. - 520.) +
               10.7 * (x1 - 520.) / (1480. - 520.);
          solver.thetal = 298.7 * (1480. - x1) / (1480. - 520.) +
                          302.4 * (x1 - 520.) / (1480. - 520.);
        } else if (x1 < 2000.) {
          qt = 10.7 * (2000. - x1) / (2000. - 1480.) +
               4.2 * (x1 - 1480.) / (2000. - 1480.);
          solver.thetal = 302.4 * (2000. - x1) / (2000. - 1480.) +
                          308.2 * (x1 - 1480.) / (2000. - 1480.);
        } else if (x1 < 3000.) {
          qt = 4.2 * (3000. - x1) / (3000. - 2000.) +
               3.0 * (x1 - 2000.) / (3000. - 2000.);
          solver.thetal = 308.2 * (3000. - x1) / (3000. - 2000.) +
                          311.85 * (x1 - 2000.) / (3000. - 2000.);
        }
        w2[1][iH2O] = qt / 1.E3;
        int err =
            root(rho / 2, rho, 1.E-4, &rho, solveDensityGivenThetaL, &solver);
        if (err) {
          msg << "### FATAL ERROR in ProblemGenerator::Bomex" << std::endl
              << "root solver doesn't converge" << std::endl
              << solveDensityGivenThetaL(rho / 2, &solver) << " "
              << solveDensityGivenThetaL(rho, &solver);
        }
        w2[1][IDN] = rho;
        w2[1][IPR] = w2[0][IPR] - 0.5 * (w2[0][IDN] + rho) * grav * solver.dz;

        // set T,P and composition
        phydro->w(IDN, k, j, i) = (w2[0][IDN] + w2[1][IDN]) / 2.;
        phydro->w(IPR, k, j, i) = (w2[0][IPR] + w2[1][IPR]) / 2.;
        phydro->w(iH2O, k, j, i) = (w2[0][iH2O] + w2[1][iH2O]) / 2.;

        x1 = pcoord->x1v(i);
        if (x1 < 700.) {
          phydro->w(IM2, k, j, i) = -8.75;
        } else if (x1 < 3000.) {
          phydro->w(IM2, k, j, i) = (-8.75) * (3000. - x1) / (3000. - 700.) +
                                    (-4.61) * (x1 - 700.) / (3000. - 700.);
        }
        phydro->w(IM1, k, j, i) = 0.;
        phydro->w(IM3, k, j, i) = 0.;

        //! \todo add noise
        //! \todo add tke
      }
    }

  // Change primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
  FreeCArray(w2);
}
