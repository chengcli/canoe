// C/C++
#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// climath
#include <climath/root.h>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/stride_iterator.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "gaussian_process.hpp"
#include "profile_inversion.hpp"

ProfileInversion::~ProfileInversion() {}

ProfileInversion::ProfileInversion(MeshBlock *pmb, ParameterInput *pin,
                                   std::string name)
    : Inversion(pmb, pin, name) {
  Application::Logger app("inversion");
  app->Log("Initializing ProfileInversion");
  char buf[80];

  // species id
  idx_ =
      Vectorize<int>(pin->GetString("inversion", name + ".variables").c_str());

  // read in prior
  for (auto m : idx_) {
    if (m == IDN) {  // change temperature
      Xstd_[IDN] = pin->GetReal("inversion", name + ".tem.std");
      Xlen_[IDN] =
          pin->GetReal("inversion", name + ".tem.corr.km") * 1.E3;  // km -> m
      app->Log(name + "::temperature std" + std::to_string(Xstd_[IDN]));
      app->Log(name + "::temperature correlation length" +
               std::to_string(Xlen_[0]));
    } else {
      Xstd_[m] = pin->GetReal("inversion", name + ".qvapor" +
                                               std::to_string(m) + ".std.gkg") /
                 1.E3;
      // km -> m
      Xlen_[m] = pin->GetReal("inversion", name + ".qvapor" +
                                               std::to_string(m) + ".corr.km") *
                 1.E3;
      snprintf(buf, sizeof(buf),
               "%s::vapor %d standard deviation = ", name.c_str(), m);
      app->Log(buf + std::to_string(Xstd_[m]));
      snprintf(buf, sizeof(buf),
               "%s::vapor %d correlation length = ", name.c_str(), m);
      app->Log(buf + std::to_string(Xlen_[m]));
    }
  }

  // power law coefficient
  chi_ = pin->GetOrAddReal("inversion", name + ".chi", 0.0);

  // Pressure sample
  plevel_ =
      Vectorize<Real>(pin->GetString("inversion", name + ".PrSample").c_str());
  int ndim = idx_.size() * plevel_.size();

  app->Log(name + "::number of input dimension = " + std::to_string(ndim));
  // app->Log(name + "::inversion pressure levels (bars) = ", plevel_);

  // add boundaries
  Real pmax = pin->GetReal("inversion", name + ".Pmax");
  Real pmin = pin->GetReal("inversion", name + ".Pmin");
  if (pmax < (plevel_.front() + 1.E-6))
    throw RuntimeError("ProfileInversion",
                       "Pmax < " + std::to_string(plevel_.front()));

  if (pmin > (plevel_.back() - 1.E-6))
    throw RuntimeError("ProfileInversion",
                       "Pmin > " + std::to_string(plevel_.back()));

  plevel_.insert(plevel_.begin(), pmax);
  plevel_.push_back(pmin);
  // app->Log(name + "::top boundary = ", pmin);
  // app->Log(name + "::bottom boundary =", pmax);

  // bar -> pa
  for (size_t n = 0; n < plevel_.size(); ++n) plevel_[n] *= 1.E5;

  // output dimension
  int nvalue = target_.size();
  // app->Log("number of output dimension = ", nvalue);

  // number of walkers
  int nwalker = pmb->block_size.nx3;
  // app->Log("walkers per block =", nwalker);
  // app->Log("total number of walkers = ", pmb->pmy_mesh->mesh_size.nx3);
  if ((nwalker < 2) && pmb->pmy_mesh->nlim > 0) {
    throw RuntimeError("ProfileInversion", "nwalker (nx3) must be at least 2");
  }

  // initialize mcmc chain
  InitializeChain(pmb->pmy_mesh->nlim + 1, nwalker, ndim, nvalue);
}

void ProfileInversion::InitializePositions() {
  int nwalker = GetWalkers();
  int ndim = GetDims();
  int nsample = samples();  // excluding boundaries

  Application::Logger app("inversion");

  // initialize random positions
  app->Log("Initializing random positions for walkers");

  unsigned int seed = time(NULL) + Globals::my_rank;
  NewCArray(init_pos_, nwalker, ndim);

  Real reference_pressure = pmy_block_->pcoord->GetReferencePressure();
  for (int n = 0; n < nwalker; ++n) {
    int ip = 0;
    for (auto m : idx_) {
      for (int i = 0; i < nsample; ++i)
        init_pos_[n][ip * nsample + i] =
            (1. * rand_r(&seed) / RAND_MAX - 0.5) * Xstd_[m] *
            pow(reference_pressure / plevel_[i + 1], chi_);
      ip++;
    }
  }
}

void ProfileInversion::UpdateHydro(Hydro *phydro, ParameterInput *pin) const {
  // read profile updates from input
  char buf[80];
  int nsample = samples();  // excluding boundaries
  Application::Logger app("inversion");
  app->Log("UpdateHydro");

  // prepare inversion model
  Real **XpSample;
  NewCArray(XpSample, 1 + NVAPOR, nsample + 2);
  std::fill(*XpSample, *XpSample + (1 + NVAPOR) * (nsample + 2), 0.);

  int is = pmy_block_->is, js = pmy_block_->js, ks = pmy_block_->ks;
  int ie = pmy_block_->ie, ke = pmy_block_->ke;

  for (auto m : idx_) {
    if (m == IDN) {  // temperature
      snprintf(buf, sizeof(buf), "%s.tema.K", name_.c_str());
    } else {  // vapors
      snprintf(buf, sizeof(buf), "%s.qvapor%da.gkg", name_.c_str(), m);
    }
    std::string str = pin->GetString("problem", buf);
    std::vector<Real> qa = Vectorize<Real>(str.c_str(), " ,");

    if (qa.size() != nsample) {
      throw ValueError("UpdateHydro", buf, qa.size(), nsample);
    }

    // g/kg -> kg/kg
    for (int i = 1; i <= nsample; ++i) {
      XpSample[m][i] = qa[i - 1] / 1.E3;
    }

    for (int k = ks; k <= ke; ++k) {
      // revise atmospheric profiles
      UpdateProfiles(phydro, XpSample, k, jl_, ju_);

      // set the revised profile to js-1
      for (int n = 0; n < NHYDRO; ++n)
        for (int i = is; i <= ie; ++i) {
          phydro->w(n, k, js - 1, i) = phydro->w(n, k, ju_, i);
        }
    }
  }

  FreeCArray(XpSample);
}

void ProfileInversion::UpdateProfiles(Hydro *phydro, Real **XpSample, int k,
                                      int jl, int ju) const {
  Application::Logger app("inversion");
  app->Log("Updating profiles at k = " + std::to_string(k));

  if (ju - jl != idx_.size()) {
    throw RuntimeError("UpdateProfiles",
                       "Number of allocations in x2 should be " +
                           std::to_string(idx_.size() + 1));
  }

  int is = pmy_block_->is, ie = pmy_block_->ie;

  int nlayer = ie - is + 1;
  int nsample = plevel_.size();

  std::vector<Real> zlev(nsample);
  Real P0 = pmy_block_->pcoord->GetReferencePressure();
  Real H0 = pmy_block_->pcoord->GetPressureScaleHeight();

  for (int i = 0; i < nsample; ++i) {
    zlev[i] = -H0 * log(plevel_[i] / P0);
  }
  // app->Log("sample levels in height = ", zlev);

  // calculate the covariance matrix of T
  std::vector<Real> stdAll(nlayer);
  std::vector<Real> stdSample(nsample);
  Real **Xp;
  NewCArray(Xp, 1 + NVAPOR, nlayer);

  // copy previous jl-1 -> jl .. ju
  for (int n = 0; n < NHYDRO; ++n)
    for (int j = jl; j <= ju; ++j)
      for (int i = is; i <= ie; ++i)
        phydro->w(n, k, j, i) = phydro->w(n, k, jl - 1, i);

  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;
  Real Rd = pthermo->GetRd();

  // calculate perturbed profiles
  app->Log("Update profiles");
  auto pcoord = pmb->pcoord;
  for (size_t n = 0; n < idx_.size(); ++n) {
    int jn = jl + n;
    int m = idx_[n];

    for (int i = is; i <= ie; ++i)
      stdAll[i - is] = Xstd_[m] * pow(exp(pcoord->x1v(i) / H0), chi_);

    for (int i = 0; i < nsample; ++i)
      stdSample[i] = Xstd_[m] * pow(exp(zlev[i] / H0), chi_);

    gp_predict(SquaredExponential, Xp[m], &pcoord->x1v(is), stdAll.data(),
               nlayer, XpSample[m], zlev.data(), stdSample.data(), nsample,
               Xlen_[m]);

    for (int i = is; i <= ie; ++i) {
      Real temp = pthermo->GetTemp(pmb, k, jn, i);

      // do not alter levels lower than zlev[0] or higher than zlev[nsample-1]
      if (pcoord->x1v(i) < zlev[0] || pcoord->x1v(i) > zlev[nsample - 1])
        continue;

      // save perturbed T profile
      if (m == IDN) {
        if (temp + Xp[m][i - is] < 0.) {
          temp = 1.;  // min 1K temperature
        } else {
          temp += Xp[m][i - is];
        }
        phydro->w(IDN, k, jn, i) = phydro->w(IPR, k, jn, i) /
                                   (Rd * temp * pthermo->RovRd(pmb, k, jn, i));
      } else {  // save perturbed compositional profiles
        phydro->w(m, k, jn, i) += Xp[m][i - is];
        phydro->w(m, k, jn, i) = std::max(phydro->w(m, k, jn, i), 0.);
        phydro->w(m, k, jn, i) = std::min(phydro->w(m, k, jn, i), 1.);
      }
    }
  }

  // copy final adjusted profiles to ju
  // set density to reflect temperature
  Real *temp1d = new Real[nlayer];

  for (int i = is; i <= ie; ++i) {
    // get initial temperature from jl-1
    temp1d[i - is] = pthermo->GetTemp(pmb, k, jl - 1, i);

    for (size_t n = 0; n < idx_.size(); ++n) {
      int jn = jl + n;
      int m = idx_[n];
      if (m == IDN) {
        // override with new temperature
        temp1d[i - is] = pthermo->GetTemp(pmb, k, jn, i);
      } else {
        phydro->w(m, k, ju, i) = phydro->w(m, k, jn, i);
      }
    }

    // reset temperature given new concentration
    phydro->w(IDN, k, ju, i) =
        phydro->w(IPR, k, ju, i) /
        (Rd * temp1d[i - is] * pthermo->RovRd(pmb, k, ju, i));
  }

  delete[] temp1d;
  FreeCArray(Xp);

  ConvectiveAdjustment(phydro, k, ju);
}

void ProfileInversion::ConvectiveAdjustment(Hydro *phydro, int k,
                                            int ju) const {
  Application::Logger app("inversion");
  app->Log("doing convective adjustment");

  int is = pmy_block_->is, ie = pmy_block_->ie;

  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;

  auto &&ac = AirParcelHelper::gather_from_primitive(pmb, k, ju, is, ie - 1);

  for (int i = is + 1; i <= ie; ++i) {
    auto &air = ac[i - is - 1];
    air.ToMoleFraction();

    Real dlnp = log(phydro->w(IPR, k, ju, i) / phydro->w(IPR, k, ju, i - 1));
    pthermo->Extrapolate(&air, dlnp, Thermodynamics::Method::NeutralStability);

    air.ToMassFraction();

    // stability
    phydro->w(IDN, k, ju, i) = std::min(air.w[IDN], phydro->w(IDN, k, ju, i));

    // saturation
    for (int n = 1; n <= NVAPOR; ++n) {
      Real rh = pthermo->RelativeHumidity(pmb, n, k, ju, i);
      if (rh > 1.) phydro->w(n, k, ju, i) /= rh;
    }
  }
}
