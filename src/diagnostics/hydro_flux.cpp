// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// canoe
#include <configure.hpp>

// diagnostics
#include "diagnostics.hpp"

HydroFlux::HydroFlux(MeshBlock *pmb)
    : Diagnostics(pmb, "hydroflux"), ncycle_(0) {
  type = "VECTORS";
  std::string varname = "v1rho,";

  for (int n = 1; n <= NVAPOR; ++n) varname += "v1q" + std::to_string(n) + ",";
  for (int n = 0; n < 3; ++n) varname += "v1v" + std::to_string(n + 1) + ",";
  varname += "v1T";

  SetName(varname);
  data.NewAthenaArray(NHYDRO, ncells1_);

  Regroup(pmb, X1DIR);

  // calculate total volume
  total_vol_.ZeroClear();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);

      for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
        total_vol_(i) += vol_(i);
      }
    }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, total_vol_.data(), total_vol_.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, mpi_comm_);
#endif
}

void HydroFlux::Progress(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  auto pthermo = Thermodynamics::GetInstance();

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle_ == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  // sum over horizontal grids weighted by volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);

      for (int i = is; i <= ie; ++i) {
        data(0, i) += vol_(i) * w(IDN, k, j, i) * w(IVX, k, j, i);
        for (int n = 1; n < IPR; ++n)
          data(n, i) += vol_(i) * w(n, k, j, i) * w(IVX, k, j, i);
        data(IPR, i) +=
            vol_(i) * pthermo->GetTemp(pmb, k, j, i) * w(IVX, k, j, i);
      }
    }

  ncycle_++;
}

void HydroFlux::Finalize(MeshBlock *pmb) {
  // take time and spatial average
  if (ncycle_ > 0) {
    // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL,
                  MPI_SUM, mpi_comm_);
#endif

    for (int n = 0; n < NHYDRO; ++n) {
      for (int i = pmb->is; i <= pmb->ie; ++i)
        data(n, i) /= ncycle_ * total_vol_(i);
    }
  }

  // clear cycle;
  ncycle_ = 0;
}
