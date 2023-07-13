// C/C++
#include <cmath>
#include <iostream>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>

// climath
#include <climath/special.h>

// utils
#include <utils/communicator.hpp>

// harp
#include "radiation.hpp"
#include "rt_solvers.hpp"

void RadiationBand::RTSolverDisort::CalBandFlux(Direction const &rayInput,
                                                Real dist, int k, int j, int il,
                                                int iu) {
  auto pcomm = Communicator::GetInstance();
  auto pmb = pmy_band_->pmy_block_;
  auto pcoord = pmb->pcoord;
  auto planet = pmb->pimpl->prad->GetPlanet();

  auto &bflxup = pmy_band_->bflxup;
  auto &bflxdn = pmy_band_->bflxdn;
  auto &bpmom = pmy_band_->bpmom;

  auto &bflags = pmy_band_->bflags_;
  auto &temf = pmy_band_->temf_;

  auto &wmin = pmy_band_->wmin_;
  auto &wmax = pmy_band_->wmax_;
  auto &spec = pmy_band_->spec_;

  auto &tau = pmy_band_->tau_;
  auto &ssa = pmy_band_->ssa_;
  auto &pmom = pmy_band_->pmom_;

  if (ds.flag.ibcnd != 0) {
    throw ValuError("RTSolverDisort::CalRadtranFlux", "ibcnd", ds.flag.ibcnd,
                    0);
  }
  pcomm->SetColor(pmb, X1DIR);

  int nblocks = pmb->pmy_mesh->mesh_size.nx1 / pmb->block_size.nx1;
  Real *bufrecv = new Real[(iu - il) * nblocks * (ds.nmom_nstr + 3)];

  if (ds.flag.planck) {
    pcomm->GatherData(temf_ + il, bufrecv, iu - il + 1);

    for (int i = 0; i < (iu - il + 1) * nblocks; ++i) {
      int m = i / (iu - il + 1);
      ds.temper[m * (iu - il) + i % (iu - il + 1)] = bufrecv[i];
    }

    std::reverse(ds.temper, ds.temper + ds.nlyr + 1);
  }
  // for (int i = 0; i <= ds.nlyr; ++i)
  //   std::cout << ds.temper[i] << std::endl;

  ds.bc.umu0 = rayInput.mu > 1.E-3 ? rayInput.mu : 1.E-3;
  ds.bc.phi0 = rayInput.phi;
  if (ds.flag.planck) {
    ds.bc.btemp = ds.temper[ds.nlyr];
    ds.bc.ttemp = ds.temper[0];
  }

  // reset flx of this column
  for (int i = il; i <= iu; ++i) bflxup(k, j, i) = bflxdn(k, j, i) = 0.;

  AthenaArray<Real> farea(iu + 1), vol(iu + 1);
  pcoord->Face1Area(k, j, il, iu, farea);
  pcoord->CellVolume(k, j, il, iu, vol);

  if (bflags & RadiationFlags::CorrelatedK) {
    // stellar source function
    if (bflags & RadiationFlags::Star)
      ds.bc.fbeam = planet->ParentInsolationFlux(wmin, wmax, dist_au);
    // planck source function
    ds.wvnmlo = wmin;
    ds.wvnmhi = wmax;
  }

  int r = pcomm->GetRank(pmb, X1DIR);

  int npmom = bpmom.GetDim4() - 1;
  int dsize = (npmom + 3) * (iu - il);
  Real *bufsend = new Real[dsize];

  // loop over bins in the band
  for (int n = 0; n < pmy_band_->GetNumBins(); ++n) {
    if (!(bflags & RadiationFlags::CorrelatedK)) {
      // stellar source function
      if (bflags & RadiationFlags::Star)
        ds.bc.fbeam =
            planet->ParentInsolationFlux(spec[n].wav1, spec[n].wav2, dist_au);
      // planck source function
      ds.wvnmlo = spec[n].wav1;
      ds.wvnmhi = spec[n].wav2;
    }

    // pack data
    packSpectralProperties(bufsend, &tau(n, il), &ssa(n, il), &pmom(n, il, 0),
                           iu - il, npmom + 1);
    pcomm->GatherData(bufsend, bufrecv, dsize);
    unpackSpectralProperties(ds.dtauc, ds.ssalb, ds.pmom, bufrecv, iu - il,
                             npmom + 1, nblocks, ds.nmom_nstr + 1);

    // absorption
    std::reverse(ds.dtauc, ds.dtauc + ds.nlyr);

    // single scatering albedo
    std::reverse(ds.ssalb, ds.ssalb + ds.nlyr);

    // Legendre coefficients
    std::reverse(ds.pmom, ds.pmom + ds.nlyr * (ds.nmom_nstr + 1));
    for (int i = 0; i < ds.nlyr; ++i)
      std::reverse(ds.pmom + i * (ds.nmom_nstr + 1),
                   ds.pmom + (i + 1) * (ds.nmom_nstr + 1));

    // run disort
    c_disort(&ds, &ds_out);

    // Counting index
    // Example, il = 0, iu = 2, ds.nlyr = 6, partition in to 3 blocks
    // face id   -> 0 - 1 - 2 - 3 - 4 - 5 - 6
    // cell id   -> | 0 | 1 | 2 | 3 | 4 | 5 |
    // disort id -> 6 - 5 - 4 - 3 - 2 - 1 - 0
    // blocks    -> ---------       *       *
    //           ->  r = 0  *       *       *
    //           ->         ---------       *
    //           ->           r = 1 *       *
    //           ->                 ---------
    //           ->                   r = 2
    // block r = 0 gets, 6 - 5 - 4
    // block r = 1 gets, 4 - 3 - 2
    // block r = 2 gets, 2 - 1 - 0
    // accumulate flux from lines
    for (int i = il; i <= iu; ++i) {
      int m = ds.nlyr - (r * (iu - il) + i - il);
      /*! \bug does not work for spherical geometry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux up
      flxup_(n, i) = ds_out.rad[m].flup;

      /*! \bug does not work for spherical geomtry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux down
      flxdn_(n, i) = ds_out.rad[m].rfldir + ds_out.rad[m].rfldn;
      bflxup(k, j, i) += spec[n].wght * flxup_(n, i);
      bflxdn(k, j, i) += spec[n].wght * flxdn_(n, i);
    }

    // spherical correction by XIZ
    // xiz 2022 flux scaling so that the heating rate is the same as the
    // plane-parallel scheme volheating scaling: first calculate flux divergence
    // from DISORT using Plane-parallel in a cell then mulitpled by the cell
    // volume divided by dx1f then solve for F using F1*S1-F2*S2 = volheating
    // the top fluxes are the still the same as the plane-paralell values
    Real volh, bflxup1 = bflxup(k, j, iu), bflxdn1 = bflxdn(k, j, iu);
    for (int i = iu - 1; i >= il; --i) {
      // upward
      volh = (bflxup1 - bflxup(k, j, i)) / pcoord->dx1f(i) * vol(i);
      bflxup1 = bflxup(k, j, i);
      bflxup(k, j, i) = (bflxup(k, j, i + 1) * farea(i + 1) - volh) / farea(i);

      // downward
      volh = (bflxdn1 - bflxdn(k, j, i)) / pcoord->dx1f(i) * vol(i);
      bflxdn1 = bflxdn(k, j, i);
      bflxdn(k, j, i) = (bflxdn(k, j, i + 1) * farea(i + 1) - volh) / farea(i);
    }
  }

  delete[] bufsend;
  delete[] bufrecv;
}

void RadiationBand::RTSolverDisort::CalBandRadiance(Direction const &rayInput,
                                                    Real dist, int k, int j,
                                                    int il, int iu) {
  throw NotImplementedError("RTSolverDisort::CalBandRadiance");
}
