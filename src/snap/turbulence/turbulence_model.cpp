// C/C++ headers
#include <algorithm>
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/reconstruct/reconstruction.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// snap
#include "turbulence_model.hpp"

// constructor, initializes data structures and parameters

TurbulenceModel::TurbulenceModel(MeshBlock *pmb, ParameterInput *pin, int nvar)
    : s(nvar, pmb->ncells3, pmb->ncells2, pmb->ncells1),
      s1(nvar, pmb->ncells3, pmb->ncells2, pmb->ncells1),
      r(nvar, pmb->ncells3, pmb->ncells2, pmb->ncells1),
      s_flux{{nvar, pmb->ncells3, pmb->ncells2, pmb->ncells1 + 1},
             {nvar, pmb->ncells3, pmb->ncells2 + 1, pmb->ncells1,
              (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated
                                 : AthenaArray<Real>::DataStatus::empty)},
             {nvar, pmb->ncells3 + 1, pmb->ncells2, pmb->ncells1,
              (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated
                                 : AthenaArray<Real>::DataStatus::empty)}},
      coarse_s_(
          nvar, pmb->ncc3, pmb->ncc2, pmb->ncc1,
          (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated
                                     : AthenaArray<Real>::DataStatus::empty)),
      coarse_r_(
          nvar, pmb->ncc3, pmb->ncc2, pmb->ncc1,
          (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated
                                     : AthenaArray<Real>::DataStatus::empty)),
      sbvar(pmb, &s, &coarse_s_, s_flux, true),
      mut(pmb->ncells3, pmb->ncells2, pmb->ncells1),
      pmy_block(pmb) {
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  Mesh *pm = pmy_block->pmy_mesh;

  Application::Logger app("snap");
  app->Log("Initialize Turbulence");

  pmb->RegisterMeshBlockData(s);

  // If user-requested time integrator is type 3S*, allocate additional memory
  // registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    s2.NewAthenaArray(nvar, nc3, nc2, nc1);
  }

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmy_block->pmr->AddToRefinement(&s, &coarse_s_);
  }

  // enroll CellCenteredBoundaryVariable object
  sbvar.bvar_index = pmb->pbval->bvars.size();
  std::cout << "bvar_index = " << sbvar.bvar_index << std::endl;
  pmb->pbval->bvars.push_back(&sbvar);
  pmb->pbval->bvars_main_int.push_back(&sbvar);

  // Allocate memory for scratch arrays
  rl_.NewAthenaArray(nvar, nc1);
  rr_.NewAthenaArray(nvar, nc1);
  rlb_.NewAthenaArray(nvar, nc1);
  x1face_area_.NewAthenaArray(nc1 + 1);
  if (pm->f2) {
    x2face_area_.NewAthenaArray(nc1);
    x2face_area_p1_.NewAthenaArray(nc1);
  }
  if (pm->f3) {
    x3face_area_.NewAthenaArray(nc1);
    x3face_area_p1_.NewAthenaArray(nc1);
  }
  cell_volume_.NewAthenaArray(nc1);
  dflx_.NewAthenaArray(nvar, nc1);

  mut.ZeroClear();
  s.ZeroClear();
  s1.ZeroClear();
  r.ZeroClear();
}

//----------------------------------------------------------------------------------------
//! \fn  void Turbulence::calculateFluxes
//  \brief Calculate passive scalar fluxes using reconstruction + weighted
//  upwinding rule

void TurbulenceModel::calculateFluxes(AthenaArray<Real> &r, const int order) {
  MeshBlock *pmb = pmy_block;

  // design decision: do not pass Hydro::flux (for mass flux) via function
  // parameters, since 1) it is unlikely that anything else would be passed, 2)
  // the current TurbulenceModel class/feature implementation is inherently
  // coupled to Hydro class 3) high-order calculation of scalar fluxes will
  // require other Hydro flux approximations (flux_fc in calculate_fluxes.cpp is
  // currently not saved persistently in Hydro class but each flux dir is temp.
  // stored in 4D scratch array scr1_nkji_)

  Hydro &hyd = *(pmb->phydro);

  AthenaArray<Real> &x1flux = s_flux[X1DIR];
  AthenaArray<Real> mass_flux;
  mass_flux.InitWithShallowSlice(hyd.flux[X1DIR], 4, IDN, 1);
  int is = pmb->is;
  int js = pmb->js;
  int ks = pmb->ks;
  int ie = pmb->ie;
  int je = pmb->je;
  int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  //--------------------------------------------------------------------------------------
  // i-direction

  // set the loop limits
  jl = js, ju = je, kl = ks, ku = ke;
  // TODO(felker): fix loop limits for fourth-order hydro
  //  if (MAGNETIC_FIELDS_ENABLED) {
  if (pmb->block_size.nx2 > 1) {
    if (pmb->block_size.nx3 == 1)  // 2D
      jl = js - 1, ju = je + 1, kl = ks, ku = ke;
    else  // 3D
      jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;
  }
  //  }

  int nvar = s.GetDim4();

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is - 1, ie + 1, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is - 1, ie + 1, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is - 1, ie + 1, r, rl_, rr_);
        for (int n = 0; n < nvar; ++n) {
#pragma omp simd
          for (int i = is; i <= ie + 1; ++i) {
            // apply floor correction
            rl_(n, i) = std::max(rl_(n, i), 0.);
            rr_(n, i) = std::max(rr_(n, i), 0.);
          }
        }
      }

      computeUpwindFlux(k, j, is, ie + 1, rl_, rr_, mass_flux, x1flux);
    }
  }

  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = s_flux[X2DIR];
    mass_flux.InitWithShallowSlice(hyd.flux[X2DIR], 4, IDN, 1);

    // set the loop limits
    il = is - 1, iu = ie + 1, kl = ks, ku = ke;
    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx3 == 1)  // 2D
      kl = ks, ku = ke;
    else  // 3D
      kl = ks - 1, ku = ke + 1;
    //    }

    for (int k = kl; k <= ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js - 1, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js - 1, il, iu, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js - 1, il, iu, r, rl_, rr_);
        for (int n = 0; n < nvar; ++n) {
#pragma omp simd
          for (int i = il; i <= iu; ++i) {
            rl_(n, i) = std::max(rl_(n, i), 0.);
            // pmb->peos->ApplyPassiveScalarFloors(rr_, n, k, j, i);
          }
        }
      }
      for (int j = js; j <= je + 1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, r, rlb_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, r, rlb_, rr_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, r, rlb_, rr_);
          for (int n = 0; n < nvar; ++n) {
#pragma omp simd
            for (int i = il; i <= iu; ++i) {
              rlb_(n, i) = std::max(rlb_(n, i), 0.);
              rr_(n, i) = std::max(rr_(n, i), 0.);
            }
          }
        }

        computeUpwindFlux(k, j, il, iu, rl_, rr_, mass_flux, x2flux);

        // swap the arrays for the next step
        rl_.SwapAthenaArray(rlb_);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = s_flux[X3DIR];
    mass_flux.InitWithShallowSlice(hyd.flux[X3DIR], 4, IDN, 1);

    // set the loop limits
    // TODO(felker): fix loop limits for fourth-order hydro
    //    if (MAGNETIC_FIELDS_ENABLED)
    il = is - 1, iu = ie + 1, jl = js - 1, ju = je + 1;

    for (int j = jl; j <= ju; ++j) {  // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks - 1, j, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks - 1, j, il, iu, r, rl_, rr_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks - 1, j, il, iu, r, rl_, rr_);
        for (int n = 0; n < nvar; ++n) {
#pragma omp simd
          for (int i = il; i <= iu; ++i) {
            rl_(n, ks - 1, j, i) = std::max(rl_(n, ks - 1, j, i), 0.);
            // pmb->peos->ApplyPassiveScalarFloors(rr_, n, k, j, i);
          }
        }
      }
      for (int k = ks; k <= ke + 1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, r, rlb_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, r, rlb_, rr_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, r, rlb_, rr_);
          for (int n = 0; n < nvar; ++n) {
#pragma omp simd
            for (int i = il; i <= iu; ++i) {
              rlb_(n, i) = std::max(rlb_(n, i), 0.);
              rr_(n, i) = std::max(rr_(n, i), 0.);
            }
          }
        }

        computeUpwindFlux(k, j, il, iu, rl_, rr_, mass_flux, x3flux);

        // swap the arrays for the next step
        rl_.SwapAthenaArray(rlb_);
      }
    }
  }

  return;
}

void TurbulenceModel::computeUpwindFlux(
    const int k, const int j, const int il,
    const int iu,                                  // CoordinateDirection dir,
    AthenaArray<Real> &rl, AthenaArray<Real> &rr,  // 2D
    AthenaArray<Real> &mass_flx,                   // 3D
    AthenaArray<Real> &flx_out) {                  // 4D
  for (int n = 0; n < flx_out.GetDim4(); n++) {
#pragma omp simd
    for (int i = il; i <= iu; i++) {
      Real fluid_flx = mass_flx(k, j, i);
      if (fluid_flx >= 0.0)
        flx_out(n, k, j, i) = fluid_flx * rl(n, i);
      else
        flx_out(n, k, j, i) = fluid_flx * rr(n, i);
    }
  }
  return;
}

void TurbulenceModel::ConservedToPrimitive(AthenaArray<Real> &s,
                                           AthenaArray<Real> const &w,
                                           AthenaArray<Real> &r,
                                           Coordinates *pco, int il, int iu,
                                           int jl, int ju, int kl, int ku) {
  int nvar = s.GetDim4();
  for (int n = 0; n < nvar; ++n)
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i) {
          r(n, k, j, i) = s(n, k, j, i) / w(IDN, k, j, i);
        }
}

void TurbulenceModel::PrimitiveToConserved(AthenaArray<Real> &r,
                                           AthenaArray<Real> const &w,
                                           AthenaArray<Real> &s,
                                           Coordinates *pco, int il, int iu,
                                           int jl, int ju, int kl, int ku) {
  int nvar = s.GetDim4();
  for (int n = 0; n < nvar; ++n)
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i) {
          s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
        }
}

void TurbulenceModel::applyBoundaryCondition(AthenaArray<Real> &r,
                                             AthenaArray<Real> &s,
                                             AthenaArray<Real> const &w,
                                             Coordinates *pco) {
  MeshBlock *pmb = pmy_block;
  int nvar = r.GetDim4();
  int bis = pmb->is - NGHOST, bie = pmb->ie + NGHOST, bjs = pmb->js,
      bje = pmb->je, bks = pmb->ks, bke = pmb->ke;

  // inner x1
  if (pmb->pbval->isPhysicalBoundary(BoundaryFace::inner_x1)) {
    for (int n = 0; n < nvar; ++n)
      for (int k = bks; k <= bke; ++k)
        for (int j = bjs; j <= bje; ++j)
          for (int i = pmb->is - NGHOST; i <= pmb->is - 1; ++i) {
            r(n, k, j, i) = r(n, k, j, pmb->is);
            s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
          }
  }

  // outer x1
  if (pmb->pbval->isPhysicalBoundary(BoundaryFace::outer_x1)) {
    for (int n = 0; n < nvar; ++n)
      for (int k = bks; k <= bke; ++k)
        for (int j = bjs; j <= bje; ++j)
          for (int i = pmb->ie + 1; i <= pmb->ie + NGHOST; ++i) {
            r(n, k, j, i) = r(n, k, j, pmb->ie);
            s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
          }
  }

  if (pmb->pmy_mesh->f2) {
    // inner x2
    if (pmb->pbval->isPhysicalBoundary(BoundaryFace::inner_x2)) {
      for (int n = 0; n < nvar; ++n)
        for (int k = bks; k <= bke; ++k)
          for (int j = pmb->js - NGHOST; j <= pmb->js - 1; ++j)
            for (int i = bis; i <= bie; ++i) {
              r(n, k, j, i) = r(n, k, pmb->js, i);
              s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
            }
    }

    // outer x2
    if (pmb->pbval->isPhysicalBoundary(BoundaryFace::outer_x2)) {
      for (int n = 0; n < nvar; ++n)
        for (int k = bks; k <= bke; ++k)
          for (int j = pmb->je + 1; j <= pmb->je + NGHOST; ++j)
            for (int i = bis; i <= bie; ++i) {
              r(n, k, j, i) = r(n, k, pmb->je, i);
              s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
            }
    }
  }

  if (pmb->pmy_mesh->f3) {
    bjs = pmb->js - NGHOST;
    bje = pmb->je + NGHOST;

    // inner x3
    if (pmb->pbval->isPhysicalBoundary(BoundaryFace::inner_x3)) {
      for (int n = 0; n < nvar; ++n)
        for (int k = pmb->ks - NGHOST; k <= pmb->js - 1; ++k)
          for (int j = bjs; j <= bje; ++j)
            for (int i = bis; i <= bie; ++i) {
              r(n, k, j, i) = r(n, pmb->ks, j, i);
              s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
            }
    }

    // outer x3
    if (pmb->pbval->isPhysicalBoundary(BoundaryFace::outer_x3)) {
      for (int n = 0; n < nvar; ++n)
        for (int k = pmb->ke + 1; k <= pmb->ke + NGHOST; ++k)
          for (int j = bjs; j <= bje; ++j)
            for (int i = bis; i <= bie; ++i) {
              r(n, k, j, i) = r(n, pmb->ke, j, i);
              s(n, k, j, i) = r(n, k, j, i) * w(IDN, k, j, i);
            }
    }
  }
}

size_t TurbulenceModel::getRestartDataSizeInBytes() {
  return s.GetSizeInBytes();
}

size_t TurbulenceModel::dumpRestartData(char *pdst) {
  std::memcpy(pdst, s.data(), s.GetSizeInBytes());
  return s.GetSizeInBytes();
}

size_t TurbulenceModel::loadRestartData(char *psrc) {
  std::memcpy(s.data(), psrc, s.GetSizeInBytes());
  // load it into the other memory register(s) too
  std::memcpy(s1.data(), psrc, s1.GetSizeInBytes());
  return s.GetSizeInBytes();
}
