// C/C++
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>
#include <athena/utils/utils.hpp>

// canoe
#include <configure.hpp>

// climath
#include <climath/core.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/extract_substring.hpp>
#include <utils/fileio.hpp>
#include <utils/ndarrays.hpp>

// opacity
#include <opacity/absorber.hpp>

// harp
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "rt_solvers.hpp"

RadiationBand::RadiationBand(std::string myname, YAML::Node const &rad)
    : NamedGroup(myname) {
  Application::Logger app("harp");
  app->Log("Initialize RadiationBand " + myname);

  if (!rad[myname]) {
    throw NotFoundError("RadiationBand", myname);
  }

  auto my = rad[myname];

  if (my["parameters"]) {
    SetRealsFrom(my["parameters"]);
  }

  pgrid_ = SpectralGridFactory::CreateFrom(my);

  wrange_ = pgrid_->ReadRangeFrom(my);

  if (my["outdir"]) {
    if (!my["outdir"].IsSequence()) {
      throw RuntimeError("RadiationBand", "outdir must be a sequence");
    }

    for (const auto &item : my["outdir"]) {
      rayOutput_.push_back(
          RadiationHelper::parse_radiation_direction(item.as<std::string>()));
    }
  }

  // set absorbers
  if (my["opacity"]) {
    if (!my["opacity"].IsSequence()) {
      throw RuntimeError("RadiationBand", "opacity must be a sequence");
    }

    auto names = my["opacity"].as<std::vector<std::string>>();
    absorbers_ = AbsorberFactory::CreateFrom(names, GetName(), rad);

    for (auto &ab : absorbers_) {
      ab->LoadOpacity(RadiationBandsFactory::GetBandId(myname));
    }
  }

  // set flags
  if (my["flags"]) {
    if (!my["flags"].IsSequence()) {
      throw RuntimeError("RadiationBand", "flags must be a sequence");
    }

    auto flag_strs = my["flags"].as<std::vector<std::string>>();
    for (auto flag : flag_strs) {
      SetFlag(RadiationHelper::parse_radiation_flags(flag));
    }
  }

  // set rt solver
  if (my["rt-solver"]) {
    psolver_ = CreateRTSolverFrom(my["rt-solver"].as<std::string>(), rad);
  } else {
    psolver_ = nullptr;
  }

  char buf[80];
  snprintf(buf, sizeof(buf), "%.2f - %.2f", wrange_.first, wrange_.second);
  app->Log("Spectral range", buf);
  app->Log("Number of spectral bins", pgrid_->spec.size());
}

RadiationBand::~RadiationBand() {
  Application::Logger app("harp");
  app->Log("Destroy RadiationBand " + GetName());
}

void RadiationBand::Resize(int nc1, int nc2, int nc3, int nstr) {
  // allocate memory for spectral properties
  tem_.resize(nc1);
  temf_.resize(nc1 + 1);

  tau_.NewAthenaArray(pgrid_->spec.size(), nc1);
  tau_.ZeroClear();

  ssa_.NewAthenaArray(pgrid_->spec.size(), nc1);
  ssa_.ZeroClear();

  pmom_.NewAthenaArray(pgrid_->spec.size(), nc1, nstr + 1);
  pmom_.ZeroClear();

  // spectral grids properties
  toa_.NewAthenaArray(pgrid_->spec.size(), rayOutput_.size(), nc3, nc2);
  toa_.ZeroClear();

  flxup_.NewAthenaArray(pgrid_->spec.size(), nc3, nc2, nc1);
  flxup_.ZeroClear();

  flxdn_.NewAthenaArray(pgrid_->spec.size(), nc3, nc2, nc1);
  flxdn_.ZeroClear();

  // band properties
  btau.NewAthenaArray(nc3, nc2, nc1);
  bssa.NewAthenaArray(nc3, nc2, nc1);
  bpmom.NewAthenaArray(nstr + 1, nc3, nc2, nc1);

  btoa.NewAthenaArray(rayOutput_.size(), nc3, nc2);
  bflxup.NewAthenaArray(nc3, nc2, nc1 + 1);
  bflxdn.NewAthenaArray(nc3, nc2, nc1 + 1);

  if (psolver_ != nullptr) {
    psolver_->Resize(nc1 - 2 * NGHOST, nstr);
  }
}

AbsorberPtr RadiationBand::GetAbsorberByName(std::string const &name) {
  for (auto &absorber : absorbers_) {
    if (absorber->GetName() == name) {
      return absorber;
    }
  }

  throw NotFoundError("RadiationBand", "Absorber " + name);

  return nullptr;
}

RadiationBand const *RadiationBand::CalBandFlux(MeshBlock const *pmb, int k,
                                                int j, int il, int iu) {
  // reset flux of this column
  for (int i = il; i <= iu; ++i) {
    bflxup(k, j, i) = 0.;
    bflxdn(k, j, i) = 0.;
  }

  psolver_->Prepare(pmb, k, j);
  psolver_->CalBandFlux(pmb, k, j, il, iu);

  return this;
}

RadiationBand const *RadiationBand::CalBandFluxColumn(MeshBlock const *pmb,
                                                      int k, int j) {
  int il = NGHOST;
  int iu = NGHOST + GetNumLayers();

  return CalBandFlux(pmb, k, j, il, iu);
}

RadiationBand const *RadiationBand::CalBandRadiance(MeshBlock const *pmb, int k,
                                                    int j) {
  // reset radiance of this column
  for (int n = 0; n < GetNumOutgoingRays(); ++n) {
    btoa(n, k, j) = 0.;
  }

  psolver_->Prepare(pmb, k, j);
  psolver_->CalBandRadiance(pmb, k, j);

  return this;
}

void RadiationBand::WriteAsciiHeader(OutputParameters const *pout) const {
  if (!TestFlag(RadiationFlags::WriteBinRadiance)) return;

  char fname[80], number[6];
  snprintf(number, sizeof(number), "%05d", pout->file_number);
  snprintf(fname, sizeof(fname), "%s.radiance.%s.txt", GetName().c_str(),
           number);
  FILE *pfile = fopen(fname, "w");

  fprintf(pfile, "# Bin Radiances of Band %s: %.3g - %.3g\n", GetName().c_str(),
          wrange_.first, wrange_.second);
  fprintf(pfile, "# Ray output size: %lu\n", rayOutput_.size());

  fprintf(pfile, "# Polar angles: ");
  for (auto &ray : rayOutput_) {
    fprintf(pfile, "%.3f", rad2deg(acos(ray.mu)));
  }
  fprintf(pfile, "\n");

  fprintf(pfile, "# Azimuthal angles: ");
  for (auto &ray : rayOutput_) {
    fprintf(pfile, "%.3f", rad2deg(ray.phi));
  }
  fprintf(pfile, "\n");

  fprintf(pfile, "#%12s%12s", "wave1", "wave2");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%12s%lu", "Radiance", j + 1);
  }
  fprintf(pfile, "%12s\n", "weight");

  fclose(pfile);
}

void RadiationBand::WriteAsciiData(OutputParameters const *pout) const {
  if (!TestFlag(RadiationFlags::WriteBinRadiance)) return;

  char fname[80], number[6];
  snprintf(number, sizeof(number), "%05d", pout->file_number);
  snprintf(fname, sizeof(fname), "%s.radiance.%s.txt", GetName().c_str(),
           number);
  FILE *pfile = fopen(fname, "w");

  for (size_t i = 0; i < pgrid_->spec.size(); ++i) {
    fprintf(pfile, "%13.3g%12.3g", pgrid_->spec[i].wav1, pgrid_->spec[i].wav2);
    for (size_t j = 0; j < rayOutput_.size(); ++j) {
      fprintf(pfile, "%12.3f", toa_(i, j));
    }
    if (TestFlag(RadiationFlags::Normalize) &&
        (wrange_.first != wrange_.second)) {
      fprintf(pfile, "%12.3g\n",
              pgrid_->spec[i].wght / (wrange_.second - wrange_.first));
    } else {
      fprintf(pfile, "%12.3g\n", pgrid_->spec[i].wght);
    }
  }

  fclose(pfile);
}

std::string RadiationBand::ToString() const {
  std::stringstream ss;
  ss << "RadiationBand: " << GetName() << std::endl;
  ss << "Absorbers: ";
  for (auto &ab : absorbers_) {
    ss << ab->GetName() << ", ";
  }
  ss << std::endl << "RT-Solver: " << psolver_->GetName();
  return ss.str();
}

std::shared_ptr<RadiationBand::RTSolver> RadiationBand::CreateRTSolverFrom(
    std::string const &rt_name, YAML::Node const &rad) {
  std::shared_ptr<RTSolver> psolver;

  if (rt_name == "Lambert") {
    psolver = std::make_shared<RTSolverLambert>(this, rad);
#ifdef RT_DISORT
  } else if (rt_name == "Disort") {
    psolver = std::make_shared<RTSolverDisort>(this, rad);
#endif  // RT_DISORT
  } else {
    throw NotFoundError("RadiationBand", rt_name);
  }

  return psolver;
}
