// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/hydro/srcterms/hydro_srcterms.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// canoe
#include "diagnostics.hpp"

LatentHeating::LatentHeating(MeshBlock *pmb)
    : Diagnostics(pmb, "Lvheating", "Latent Heating") {
  type = "SCALARS";
  grid = "--C";
  units = "W/m^3";
  data.NewAthenaArray(1, 1, 1, ncells1_);
}

LatentHeating::~LatentHeating() { data.DeleteAthenaArray(); }

void LatentHeating::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  int iH2O = 1, iH2Oc = 3, iH2Op = 5;

  AthenaArray<Real> &x1flux = phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux = phydro->flux[X2DIR];
  AthenaArray<Real> &x3flux = phydro->flux[X3DIR];

  if (ncycle == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        // mass flux divergence of water cloud, and droplets
        Real div = 0.;
        // surface area
        Real x1area_l, x1area_r, x2area_l, x2area_r, x3area_l, x3area_r;
        // calculate area of each interface
        x1area_l = pcoord->GetFace1Area(k, j, i);
        x1area_r = pcoord->GetFace1Area(k, j, i + 1);
        x2area_l = pcoord->GetFace2Area(k, j, i);
        x2area_r = pcoord->GetFace2Area(k, j + 1, i);
        x3area_l = pcoord->GetFace3Area(k, j, i);
        x3area_r = pcoord->GetFace3Area(k + 1, j, i);
        // calculate local mass flux divergence of water [kg/(m^3 s)]
        div = (x1area_r * x1flux(iH2O, k, j, i + 1) -
               x1area_l * x1flux(iH2O, k, j, i)) +
              (x2area_r * x2flux(iH2O, k, j + 1, i) -
               x2area_l * x2flux(iH2O, k, j, i)) +
              (x3area_r * x3flux(iH2O, k + 1, j, i) -
               x3area_l * x3flux(iH2O, k, j, i));
        // data(n,k,j,i) = -div/pcoord->GetCellVolume(k,j,i);
        div /= -vol_(i);
        // calculate vapor, cloud, and precip changes, dq/dt [kg/(m^3 s)]
        Real dqdt = (w(iH2O, k, j, i) * w(IDN, k, j, i) -
                     pmb->ruser_meshblock_data[4](iH2O, k, j, i)) /
                    pmb->pmy_mesh->dt;
        // local tendency due to chemistry (evaporation and autoconversion)
        // data(n,k,j,i) = dqdt-data(n,k,j,i);
        Real Lv =
            2.2E6;  // pthermo->GetLatent(iH2O,pthermo->Temp(w.at(k,j,i)));
        data(i) -= Lv * (dqdt - div) * vol_(i);
        pmb->ruser_meshblock_data[4](iH2O, k, j, i) =
            w(iH2O, k, j, i) * w(IDN, k, j, i);
      }
    }
  ncycle++;
}

void LatentHeating::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *total_vol = new Real[ncells1_];
  std::fill(total_vol, total_vol + ncells1_, 0.);
  // calculate total volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) total_vol[i] += vol_(i);
    }

  // take time and spatial average
  if (ncycle > 0) {
    // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, total_vol, ncells1_, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    for (int i = is; i <= ie; ++i) {
      data(i) /= ncycle * total_vol[i];
    }
  }

  // clear cycle;
  delete[] total_vol;
  ncycle = 0;
}
