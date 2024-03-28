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

ProfileInversion::ProfileInversion(MeshBlock *pmb, YAML::Node const &node)
    : Inversion(pmb, "profile") {
  Application::Logger app("inversion");
  app->Log("Initializing ProfileInversion");
  char buf[80];

  // species id
  auto species = node["variables"].as<std::vector<std::string>>();
  SetSpeciesIndex(species);

  // Pressure sample
  plevel_ = node["sample-pressures"].as<std::vector<Real>>();

  // add boundaries
  Real pmax = node["bottom-pressure"].as<Real>();
  Real pmin = node["top-pressure"].as<Real>();

  if (pmax < (plevel_.front() + 1.E-6))
    throw RuntimeError("ProfileInversion",
                       "Pmax < " + std::to_string(plevel_.front()));

  if (pmin > (plevel_.back() - 1.E-6))
    throw RuntimeError("ProfileInversion",
                       "Pmin > " + std::to_string(plevel_.back()));

  plevel_.insert(plevel_.begin(), pmax);
  plevel_.push_back(pmin);

  // bar -> pa
  for (size_t n = 0; n < plevel_.size(); ++n) {
    plevel_[n] *= 1.E5;
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
      snprintf(buf, sizeof(buf), "%s.tema.K", GetName().c_str());
    } else {  // vapors
      snprintf(buf, sizeof(buf), "%s.qvapor%da.gkg", GetName().c_str(), m);
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

void ProfileInversion::UpdateModel(std::vector<Real> const &par, int k) const {
  app->Log("Update profiles");
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;

  int is = pmy_block_->is, ie = pmy_block_->ie;
  int nlayer = ie - is + 1;
  int nsample = plevel_.size();

  std::vector<Real> zlev(nsample);
  Real P0 = pmy_block_->pcoord->GetReferencePressure();
  Real H0 = pmy_block_->pcoord->GetPressureScaleHeight();

  for (int i = 0; i < nsample; ++i) {
    zlev[i] = -H0 * log(plevel_[i] / P0);
  }

  // calculate the covariance matrix of T
  std::vector<Real> stdAll(nlayer);
  std::vector<Real> stdSample(nsample);

  AirColumn ac(nlayer);

  Real **Xp;
  NewCArray(Xp, 1 + NVAPOR, nlayer);

  // copy previous jl-1 -> jl .. ju
  for (int n = 0; n < NHYDRO; ++n)
    for (int j = jl; j <= ju; ++j)
      for (int i = is; i <= ie; ++i)
        phydro->w(n, k, j, i) = phydro->w(n, k, jl - 1, i);

  Real Rd = pthermo->GetRd();

  Real chi = GetPar<Real>("chi");

  // calculate perturbed profiles
  auto pcoord = pmb->pcoord;
  for (size_t n = 0; n < GetMySpeciesIndices().size(); ++n) {
    int jn = jl + n;
    int m = mySpeciesId(n);

    for (int i = is; i <= ie; ++i)
      stdAll[i - is] = Xstd_[n] * pow(exp(pcoord->x1v(i) / H0), chi);

    for (int i = 0; i < nsample; ++i)
      stdSample[i] = Xstd_[n] * pow(exp(zlev[i] / H0), chi);

    gp_predict(SquaredExponential, Xp[m], &pcoord->x1v(is), stdAll.data(),
               nlayer, XpSample[m], zlev.data(), stdSample.data(), nsample,
               Xlen_[n]);

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
      int m = mySpeciesId(n);
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

  enforceStability(phydro, k);
}

void ProfileInversion::enforceStability(AirColumn &ac, int k) const {
  Application::Logger app("inversion");
  app->Log("Enforce stability");

  int is = pmy_block_->is, ie = pmy_block_->ie;

  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;

  for (int i = 1; i < ac.size(); ++i) {
    auto air = ac[i - 1];
    air.ToMoleFraction();

    Real dlnp = ac[i].w[IPR] / ac[i - 1].w[IPR];
    pthermo->Extrapolate(&air, dlnp, Thermodynamics::Method::NeutralStability);

    air.ToMassFraction();

    // stability
    ac[i].w[IDN] = std::min(air.w[IDN], ac[i].w[IDN]);

    // saturation
    for (int n = 1; n <= NVAPOR; ++n) {
      Real rh = pthermo->RelativeHumidity(pmb, n, k, ju_, i);
      if (rh > 1.) ac[i].w[n] /= rh;
    }
  }
}
