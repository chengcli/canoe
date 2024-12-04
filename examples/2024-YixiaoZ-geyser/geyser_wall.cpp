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

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>

Real H2Oratio, CO2ratio, grav;
int iH2O, iH2Oc, iCO2, iCO2c;

Real wall1_corner_x2;
Real wall1_corner_x1;
Real wall2_corner_x2;
Real wall2_corner_x1;

Real Ptriple1, Ttriple1;
Real Rd, gammad;
Real x1min, x1max, x2min, x2max;

Real massflux_H2ratio, massflux_CO2ratio;

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
        // add CO2
        drhoCO2 = drhoH2O * massflux_CO2ratio;
        u(iCO2, k, j, is) += drhoCO2;
        u(IEN, k, j, is) += drhoCO2 * (Rd / (gammad - 1.)) *
                            pthermo->GetCvRatio(iCO2) * Ttriple1;
      }
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  BottomInjection(pmb, time, dt, w, r, bcc, u, s);
  // WallInteraction(pmb, time, dt, w, r, bcc, u, s);
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

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  // construct atmosphere from bottom up
  std::vector<Real> yfrac(IVX, 0.);
  yfrac[iH2O] = H2Oratio;
  // yfrac[iCO2] = CO2ratio;
  yfrac[0] = 1. - H2Oratio;

  auto &w = phydro->w;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pthermo->SetMassFractions<Real>(yfrac.data());
      pthermo->EquilibrateTP(273., 1.);

      // half a grid to cell center
      pthermo->Extrapolate_inplace(pcoord->dx1f(is) / 2., "isothermal", grav);

      for (int i = is; i <= ie; ++i) {
        pthermo->GetPrimitive(w.at(k, j, i));

        // set all clouds to zero
        for (int n = 1 + NVAPOR; n < IVX; ++n) w(n, k, j, i) = 0.;

        // move to the next cell
        pthermo->Extrapolate_inplace(pcoord->dx1f(i), "isothermal", grav);
      }
    }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

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
