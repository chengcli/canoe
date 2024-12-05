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

Real H2Oratio, CO2ratio, grav;
int iH2O, iH2Oc, iCO2, iCO2c;

Real wall1_corner_x2;
Real wall1_corner_x1;
Real wall2_corner_x2;
Real wall2_corner_x1;
Real sigtanh;

Real Ptriple1, Ttriple1;
Real Rd, gammad;
Real x1min, x1max, x2min, x2max;

Real massflux_H2ratio, massflux_CO2ratio;
Real Tm, Ts;

inline double SatVaporPresIdeal(double t, double p, double beta, double gamma) {
  return p * exp((1. - 1. / t) * beta - gamma * log(t));
}

double sat_vapor_p_H2O(double T) {
  double betal = 22.46, gammal = 0, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

void reflecting_x2_left(MeshBlock *pmb, Coordinates *pco,
                        AthenaArray<Real> &prim, FaceField &b, Real time,
                        Real dt, int il, int iu, int jl, int ju, int kl, int ku,
                        int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n = 0; n < NHYDRO; ++n) {
    if (n == IVY) {
      for (int k = kl; k <= ku; ++k) {
        for (int j = 1; j <= ngh; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim(IVY, k, jl - j, i) =
                -prim(IVY, k, jl + j - 1, i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k = kl; k <= ku; ++k) {
        for (int j = 1; j <= ngh; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim(n, k, jl - j, i) = prim(n, k, jl + j - 1, i);
          }
        }
      }
    }
  }
}

void reflecting_x2_right(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim, FaceField &b, Real time,
                         Real dt, int il, int iu, int jl, int ju, int kl,
                         int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n = 0; n < NHYDRO; ++n) {
    if (n == (IVY)) {
      for (int k = kl; k <= ku; ++k) {
        for (int j = 1; j <= ngh; ++j) {
#pragma omp simd
          for (int i = il; i <= iu; ++i) {
            prim(IVY, k, ju + j, i) =
                -prim(IVY, k, ju - j + 1, i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k = kl; k <= ku; ++k) {
        for (int j = 1; j <= ngh; ++j) {
#pragma omp simd
          for (int i = il; i <= iu; ++i) {
            prim(n, k, ju + j, i) = prim(n, k, ju - j + 1, i);
          }
        }
      }
    }
  }
}

void reflecting_x1_left(MeshBlock *pmb, Coordinates *pco,
                        AthenaArray<Real> &prim, FaceField &b, Real time,
                        Real dt, int il, int iu, int jl, int ju, int kl, int ku,
                        int ngh) {
  for (int n = 0; n < NHYDRO; ++n) {
    if (n == IVX) {
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = 1; i <= ngh; ++i) {
            prim(n, k, j, il - i) = -prim(n, k, j, il + i - 1);
          }
        }
      }
    } else {
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = 1; i <= ngh; ++i) {
            prim(n, k, j, il - i) = prim(n, k, j, il + i - 1);
          }
        }
      }
    }
  }
}

void reflecting_x1_right(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim, FaceField &b, Real time,
                         Real dt, int il, int iu, int jl, int ju, int kl,
                         int ku, int ngh) {
  for (int n = 0; n < NHYDRO; ++n) {
    if (n == (IVX)) {
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
#pragma omp simd
          for (int i = 1; i <= ngh; ++i) {
            prim(n, k, j, iu + i) = -prim(n, k, j, iu - i + 1);
          }
        }
      }
    } else {
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
#pragma omp simd
          for (int i = 1; i <= ngh; ++i) {
            prim(n, k, j, iu + i) = prim(n, k, j, iu - i + 1);
          }
        }
      }
    }
  }
}

bool fclose(Real x, Real x0) { return std::abs(x - x0) < 1.e-6; }

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

void WallInteraction(MeshBlock *pmb, Real const time, Real const dt,
                     AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                     AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                     AthenaArray<Real> &s) {
  int js = pmb->js;
  int je = pmb->je;
  int is = pmb->is;
  int ie = pmb->ie;
  int jw;  // index of j at the wall

  auto pthermo = Thermodynamics::GetInstance();

  Real p_H2O, drhoH2O, drhoH2, drhoCO2;
  Real Tw, Pw, Ta, z, csw, csa, KE;
  Real drag_coef = 2e-3;
  Real dv1;

  Real x2f_left, x2f_right, x1f_left, x1f_right, x1f_center, x2f_center;
  Real tanhweight;

  // remove water vapor
  for (int jj = 0; jj <= 1; ++jj) {
    jw = (jj == 0) ? js : je;
    x2f_left = pmb->pcoord->x2f(jw);
    x2f_right = pmb->pcoord->x2f(jw + 1);
    x2f_center = (x2f_left + x2f_right) / 2;
    x1f_left = pmb->pcoord->x1f(is);
    x1f_right = pmb->pcoord->x1f(ie + 1);
    x1f_center = (x1f_left + x1f_right) / 2;

    if ((x2f_center < wall1_corner_x2) || (x2f_center > wall2_corner_x2)) {
      if (x1f_center < wall1_corner_x1) {
        continue;
      }

      if (x1f_left - 2 * pmb->pcoord->dx1f(is) < wall1_corner_x1) {
        for (int k = pmb->ks; k <= pmb->ke; ++k)
          for (int j = pmb->js; j <= pmb->je; ++j) {
            Ta = pthermo->GetTemp(w.at(k, j, is));
            p_H2O = (pmb->phydro->w(IDN, k, j, is) *
                     pmb->phydro->w(iH2O, k, j, is) * Rd * Ta *
                     pthermo->GetInvMuRatio(iH2O));
            Pw = sat_vapor_p_H2O(Ts);
            csw = sqrt(2 * M_PI * Rd * Ts * pthermo->GetInvMuRatio(iH2O));
            csa = sqrt(2 * M_PI * Rd * Ta * pthermo->GetInvMuRatio(iH2O));
            // drhoH2O = dt * (Pw/csw - p_H2O/csa) / pmb->pcoord->dx1f(is);
            drhoH2O = dt * (-p_H2O / csa) / pmb->pcoord->dx1f(is);
            u(iH2O, k, j, is) += drhoH2O;
            // std::cout << "x1min" << x1f_left << "x2min" << x2f_left << "Add
            // upper wall" << std::endl;

            if (drhoH2O < 0) {
              KE = 0.5f * (pmb->phydro->w(IVX, k, j, is) *
                               pmb->phydro->w(IVX, k, j, is) +
                           pmb->phydro->w(IVY, k, j, is) *
                               pmb->phydro->w(IVY, k, j, is) +
                           pmb->phydro->w(IVZ, k, j, is) *
                               pmb->phydro->w(IVZ, k, j, is));
              u(IEN, k, j, is) +=
                  drhoH2O *
                  (KE + (Rd / (gammad - 1.)) * pthermo->GetCvRatio(iH2O) * Ta);
              u(IVZ, k, j, is) += drhoH2O * pmb->phydro->w(IVZ, k, j, is);
              u(IVY, k, j, is) += drhoH2O * pmb->phydro->w(IVY, k, j, is);
              u(IVX, k, j, is) += drhoH2O * pmb->phydro->w(IVX, k, j, is);
            } else {
              u(IEN, k, j, is) += drhoH2O * ((Rd / (gammad - 1.)) *
                                             pthermo->GetCvRatio(iH2O) * Tw);
            }
          }
      }
      continue;
    }

    if (x1f_center > wall1_corner_x1) {
      continue;
    }
    if ((x2f_left - wall1_corner_x2 < pmb->pcoord->dx2f(js)) ||
        (wall2_corner_x2 - x2f_right < pmb->pcoord->dx2f(je))) {
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          // if (pmb->pcoord->x1f(i) < 0.2 * wall2_cornerx1) {
          //   continue;
          // }
          tanhweight =
              (1. + tanh((pmb->pcoord->x1f(i) - 5. * sigtanh) / sigtanh)) / 2.0;

          dv1 = -(dt * drag_coef * pmb->phydro->w(IDN, k, jw, i) *
                  pmb->phydro->w(IVX, k, jw, i) *
                  pmb->phydro->w(IVX, k, jw, i) / pmb->pcoord->dx2f(jw));
          u(IVX, k, jw, i) += dv1;  // add drag
          u(IEN, k, jw, i) +=
              dv1 * pmb->phydro->u(IVX, k, jw, i);  // subduct energy

          Ta = pthermo->GetTemp(w.at(k, jw, i));

          p_H2O =
              (pmb->phydro->w(IDN, k, jw, i) * pmb->phydro->w(iH2O, k, jw, i) *
               Rd * Ta * pthermo->GetInvMuRatio(iH2O));
          z = pmb->pcoord->x1f(i);
          Tw = Tm * pow(Ts / Tm, (z - x1min) / (wall1_corner_x1 - x1min));
          Pw = sat_vapor_p_H2O(Tw);

          csw = sqrt(2 * M_PI * Rd * Tw * pthermo->GetInvMuRatio(iH2O));
          csa = sqrt(2 * M_PI * Rd * Ta * pthermo->GetInvMuRatio(iH2O));

          drhoH2O = dt * (Pw / csw - p_H2O / csa) / pmb->pcoord->dx2f(jw);
          drhoH2O *= tanhweight;

          u(iH2O, k, jw, i) += drhoH2O;

          if (drhoH2O < 0) {
            KE =
                0.5f *
                (pmb->phydro->w(IVX, k, jw, i) * pmb->phydro->w(IVX, k, jw, i) +
                 pmb->phydro->w(IVY, k, jw, i) * pmb->phydro->w(IVY, k, jw, i) +
                 pmb->phydro->w(IVZ, k, jw, i) * pmb->phydro->w(IVZ, k, jw, i));
            u(IEN, k, jw, i) +=
                drhoH2O *
                (KE + (Rd / (gammad - 1.)) * pthermo->GetCvRatio(iH2O) * Ta);
            u(IVZ, k, jw, i) += drhoH2O * pmb->phydro->w(IVZ, k, jw, i);
            u(IVY, k, jw, i) += drhoH2O * pmb->phydro->w(IVY, k, jw, i);
            u(IVX, k, jw, i) += drhoH2O * pmb->phydro->w(IVX, k, jw, i);
          } else {
            u(IEN, k, jw, i) += drhoH2O * ((Rd / (gammad - 1.)) *
                                           pthermo->GetCvRatio(iH2O) * Tw);
          }
        }
    }
  }
}

void BottomInjection(MeshBlock *pmb, Real const time, Real const dt,
                     AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                     AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                     AthenaArray<Real> &s) {
  int is = pmb->is;
  int ie = pmb->ie;

  auto pthermo = Thermodynamics::GetInstance();

  Real p, drhoH2O, drhoH2, drhoCO2;

  Real x1s = pmb->pcoord->x1f(is);

  if (x1s < x1min + pmb->pcoord->dx1f(is)) {
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        // std::cout << pmb->pcoord->x2v(j) << std::endl;
        //  inject at the center of the bottom boundary
        if ((pmb->pcoord->x2v(j) < wall1_corner_x2) ||
            pmb->pcoord->x2v(j) > wall2_corner_x2) {
          continue;
        }

        p = pmb->phydro->w(IPR, k, j, is);

        // add water vapor
        drhoH2O =
            (dt * (Ptriple1 - p) /
             sqrt(2 * M_PI * Rd * Ttriple1 * pthermo->GetInvMuRatio(iH2O))) /
            pmb->pcoord->dx1f(is);
        u(iH2O, k, j, is) += drhoH2O;
        u(IEN, k, j, is) += drhoH2O * (Rd / (gammad - 1.)) *
                            pthermo->GetCvRatio(iH2O) * Ttriple1;
        // add dry air (H2)
        drhoH2 = drhoH2O * massflux_H2ratio;
        u(IDN, k, j, is) += drhoH2;
        u(IEN, k, j, is) += drhoH2 * (Rd / (gammad - 1.)) * Ttriple1;

        /* add CO2
        drhoCO2 = drhoH2O * massflux_CO2ratio;
        u(iCO2, k, j, is) += drhoCO2;
        u(IEN, k, j, is) += drhoCO2 * (Rd / (gammad - 1.)) *
                            pthermo->GetCvRatio(iCO2) * Ttriple1;*/
      }
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  BottomInjection(pmb, time, dt, w, r, bcc, u, s);
  WallInteraction(pmb, time, dt, w, r, bcc, u, s);
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  H2Oratio = pin->GetReal("initialcondition", "H2Oratio");
  CO2ratio = pin->GetReal("initialcondition", "CO2ratio");

  grav = -pin->GetReal("hydro", "grav_acc1");

  // index
  iH2O = pthermo->SpeciesIndex("H2O");
  iH2Oc = pthermo->SpeciesIndex("H2O(s)");
  // iCO2 = pthermo->SpeciesIndex("CO2");
  // iCO2c = pthermo->SpeciesIndex("CO2(s)");

  wall1_corner_x1 = pin->GetReal("problem", "wall1_corner_x1");
  wall1_corner_x2 = pin->GetReal("problem", "wall1_corner_x2");
  wall2_corner_x1 = pin->GetReal("problem", "wall2_corner_x1");
  wall2_corner_x2 = pin->GetReal("problem", "wall2_corner_x2");
  sigtanh = pin->GetReal("problem", "sigtanh");

  Ptriple1 = pin->GetReal("problem", "Ptriple1");
  Ttriple1 = pin->GetReal("problem", "Ttriple1");
  Rd = pin->GetReal("problem", "Rd");
  gammad = pin->GetReal("hydro", "gamma");
  x1min = pin->GetReal("mesh", "x1min");
  x1max = pin->GetReal("mesh", "x1max");
  x2min = pin->GetReal("mesh", "x2min");
  x2max = pin->GetReal("mesh", "x2max");

  massflux_H2ratio = pin->GetReal("problem", "massflux_H2ratio");
  massflux_CO2ratio = pin->GetReal("problem", "massflux_CO2ratio");

  Tm = pin->GetReal("problem", "Tm");
  Ts = pin->GetReal("problem", "Ts");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  // construct 1d atmosphere from bottom up
  std::vector<Real> yfrac(IVX, 0.);
  yfrac[iH2O] = H2Oratio;
  // yfrac[iCO2] = CO2ratio;
  yfrac[0] = 1. - H2Oratio;

  int nx1 = pmy_mesh->mesh_size.nx1;
  Real dz = (x1max - x1min) / (nx1 - 1);
  std::cout << "nx1 = " << nx1 << std::endl;

  AthenaArray<Real> w1, z1;
  w1.NewAthenaArray(NHYDRO, nx1);

  z1.NewAthenaArray(nx1);
  z1(0) = x1min + dz / 2.;
  for (int i = 1; i < nx1; ++i) z1(i) = z1(i - 1) + dz;

  pthermo->SetMassFractions<Real>(yfrac.data());
  pthermo->EquilibrateTP(100., 1.);

  // half a grid to cell center
  pthermo->Extrapolate_inplace(dz / 2., "isothermal", grav);

  for (int i = 0; i < nx1; ++i) {
    pthermo->GetPrimitive(w1.at(i));

    // set all clouds to zero
    for (int n = 1 + NVAPOR; n < IVX; ++n) w1(n, i) = 0.;

    // move to the next cell
    pthermo->Extrapolate_inplace(dz, "isothermal", grav);
  }

  // populate to 3D mesh
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n) {
          phydro->w(n, k, j, i) =
              interp1(pcoord->x1v(i), w1.data() + n * nx1, z1.data(), nx1);
        }
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

  // add walls
  Real x1minblock = block_size.x1min;
  Real x1maxblock = block_size.x1max;
  Real x1c = (x1minblock + x1maxblock) / 2;
  Real x2minblock = block_size.x2min;
  Real x2maxblock = block_size.x2max;
  Real x2c = (x2minblock + x2maxblock) / 2;

  if (fclose(x2minblock, wall1_corner_x2) && (x1c < wall1_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::inner_x2] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::inner_x2] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::inner_x2] = true;
    std::cout << "Boundary left enrolled" << std::endl;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::inner_x2,
                                         reflecting_x2_left);
  }

  if (fclose(x2maxblock, wall1_corner_x2) && (x1c < wall1_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::outer_x2] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::outer_x2] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::outer_x2] = true;
    std::cout << "Boundary right enrolled" << std::endl;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::outer_x2,
                                         reflecting_x2_right);
  }

  if ((x2c < wall1_corner_x2) && fclose(x1maxblock, wall1_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::outer_x1] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::outer_x1] = true;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::outer_x1,
                                         reflecting_x1_right);
  }

  if ((x2c < wall1_corner_x2) && fclose(x1minblock, wall1_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::inner_x1] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::inner_x1] = true;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::inner_x1,
                                         reflecting_x1_left);
  }

  // add a block
  if (fclose(x2minblock, wall2_corner_x2) && (x1c < wall2_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::inner_x2] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::inner_x2] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::inner_x2] = true;
    std::cout << "Boundary left enrolled" << std::endl;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::inner_x2,
                                         reflecting_x2_left);
  }

  if (fclose(x2maxblock, wall2_corner_x2) && (x1c < wall2_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::outer_x2] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::outer_x2] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::outer_x2] = true;
    std::cout << "Boundary right enrolled" << std::endl;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::outer_x2,
                                         reflecting_x2_right);
  }

  if ((x2c > wall2_corner_x2) && fclose(x1maxblock, wall2_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::outer_x1] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::outer_x1] = true;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::outer_x1,
                                         reflecting_x1_right);
  }

  if ((x2c > wall2_corner_x2) && fclose(x1minblock, wall2_corner_x1)) {
    pmy_mesh->mesh_bcs[BoundaryFace::inner_x1] = BoundaryFlag::user;
    pbval->block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::user;
    pbval->apply_bndry_fn_[BoundaryFace::inner_x1] = true;
    pmy_mesh->EnrollUserBoundaryFunction(BoundaryFace::inner_x1,
                                         reflecting_x1_left);
  }
}
