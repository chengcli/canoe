// C/C++
#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

// cantera
#include <cantera/base/yaml.h>

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
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "gaussian_process.hpp"
#include "inversion.hpp"

ProfileInversion::ProfileInversion(YAML::Node const &node)
    : Inversion("profile") {
  Application::Logger app("inversion");
  app->Log("Initializing ProfileInversion");

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

void ProfileInversion::UpdateModel(MeshBlock *pmb, std::vector<Real> const &par,
                                   int k) const {
  Application::Logger app("inversion");
  app->Log("Update Model");

  auto pthermo = Thermodynamics::GetInstance();

  int nlayer = pmb->ie - pmb->is + 1;
  int nvar = GetMySpeciesIndices().size();
  int nsample = plevel_.size();

  std::vector<Real> zlev(nsample);
  Real P0 = pmb->pcoord->GetReferencePressure();
  Real H0 = pmb->pcoord->GetPressureScaleHeight();

  for (int i = 0; i < nsample; ++i) {
    zlev[i] = -H0 * log(plevel_[i] / P0);
  }

  // calculate the covariance matrix of T
  std::vector<Real> stdAll(nlayer);
  std::vector<Real> stdSample(nsample);

  Real const **XpSample;
  Real **Xp;

  XpSample = new Real const *[nvar];
  for (size_t n = 0; n < GetMySpeciesIndices().size(); ++n) {
    XpSample[n] = &par[nvar * nsample];
  }

  NewCArray(Xp, nvar, nlayer);

  Real Rd = pthermo->GetRd();
  Real chi = GetPar<Real>("chi");
  int is = pmb->is, ie = pmb->ie;

  // calculate perturbed profiles
  auto pcoord = pmb->pcoord;
  auto ac = AirParcelHelper::gather_from_primitive(pmb, k, jl_, is, ie);

  for (size_t n = 0; n < nvar; ++n) {
    int m = mySpeciesId(n);

    for (int i = is; i <= ie; ++i)
      stdAll[i - is] = Xstd_[n] * pow(exp(pcoord->x1v(i) / H0), chi);

    for (int i = 0; i < nsample; ++i)
      stdSample[i] = Xstd_[n] * pow(exp(zlev[i] / H0), chi);

    gp_predict(SquaredExponential, Xp[n], &pcoord->x1v(is), stdAll.data(),
               nlayer, XpSample[n], zlev.data(), stdSample.data(), nsample,
               Xlen_[n]);

    for (int i = 0; i < nlayer; ++i) {
      ac[i].ToMoleFraction();

      // do not alter levels lower than zlev[0] or higher than zlev[nsample-1]
      if (pcoord->x1v(is + i) < zlev[0] ||
          pcoord->x1v(is + i) > zlev[nsample - 1])
        continue;

      // save perturbed T profile
      if (m == IDN) {
        if (ac[i].w[IDN] + Xp[n][i] < 0.) {
          ac[i].w[IDN] = 1.;  // min 1K temperature
        } else {
          ac[i].w[IDN] += Xp[n][i];
        }
      } else {  // save perturbed compositional profiles
        ac[i].w[m] += Xp[n][i - is];
        ac[i].w[m] = std::max(ac[i].w[m], 0.);
        ac[i].w[m] = std::min(ac[i].w[m], 1.);
      }
    }

    AirParcelHelper::distribute_to_primitive(pmb, k, jl_ + n, ac);
  }

  enforceStability(ac);
  AirParcelHelper::distribute_to_primitive(pmb, k, ju_, ac);

  delete[] XpSample;
  FreeCArray(Xp);
}

void ProfileInversion::enforceStability(AirColumn &ac) const {
  Application::Logger app("inversion");
  app->Log("Enforce stability");

  auto pthermo = Thermodynamics::GetInstance();

  for (int i = 1; i < ac.size(); ++i) {
    auto air = ac[i - 1];
    air.ToMoleFraction();

    Real dlnp = ac[i].w[IPR] / ac[i - 1].w[IPR];
    pthermo->Extrapolate(&air, dlnp, "neutral");

    air.ToMassFraction();

    // stability
    ac[i].w[IDN] = std::min(air.w[IDN], ac[i].w[IDN]);

    // saturation
    for (int n = 1; n <= NVAPOR; ++n) {
      Real rh = get_relative_humidity(ac[i], n);
      if (rh > 1.) ac[i].w[n] /= rh;
    }
  }
}
