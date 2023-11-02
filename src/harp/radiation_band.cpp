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
#include <utils/parameter_map.hpp>
#include <utils/vectorize.hpp>

// opacity
#include <opacity/absorber.hpp>

// harp
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "radiation_utils.hpp"  // readRadiationDirections
#include "rt_solvers.hpp"
#include "spectral_grid.hpp"

RadiationBand::RadiationBand(YAML::Node &rad, std::string myname, int nc1,
                             int nc2, int nc3)
    : NamedGroup(myname) {
  Application::Logger app("harp");
  app->Log("Initialize RadiationBand " + myname);

  if (!rad[myname]) {
    throw NotFoundError("RadiationBand", myname);
  }

  auto &my = rad[myname];

  if (my["parameters"]) {
    SetRealsFrom(my["parameters"]);
  }

  // number of phase function moments
  int npmom = my["phase-function-moments"] ? my["moments"].as<int>() : 1;

  wrange_ = pgrid_->ReadRangeFrom(my);

  pgrid_ = SpectralGridFactory::CreateFrom(my);

  // outgoing radiation direction (mu,phi) in degree
  if (pin->DoesParameterExist("radiation", GetName() + ".outdir")) {
    auto str = pin->GetString("radiation", GetName() + ".outdir");
    read_radiation_directions(&rayOutput_, str);
  } else if (pin->DoesParameterExist("radiation", "outdir")) {
    auto str = pin->GetString("radiation", "outdir");
    read_radiation_directions(&rayOutput_, str);
  }

  // allocate memory for spectral properties
  tem_.NewAthenaArray(nc1);
  temf_.NewAthenaArray(nc1 + 1);

  tau_.NewAthenaArray(pgrid_->GetSize(), nc1);
  tau_.ZeroClear();

  ssa_.NewAthenaArray(pgrid_->GetSize(), nc1);
  ssa_.ZeroClear();

  toa_.NewAthenaArray(pgrid_->GetSize(), rayOutput_.size());
  toa_.ZeroClear();

  pmom_.NewAthenaArray(pgrid_->GetSize(), nc1, npmom + 1);
  pmom_.ZeroClear();

  flxup_.NewAthenaArray(pgrid_->GetSize(), nc1);
  flxup_.ZeroClear();

  flxdn_.NewAthenaArray(pgrid_->GetSize(), nc1);
  flxdn_.ZeroClear();

  // band properties
  btau.NewAthenaArray(nc3, nc2, nc1);
  bssa.NewAthenaArray(nc3, nc2, nc1);
  bpmom.NewAthenaArray(npmom + 1, nc3, nc2, nc1);

  //! \note btoa, bflxup, bflxdn are shallow slices to Radiation variables

  // set absorbers
  if (my["opacity"]) {
    auto absorber_names = my["opacity"].as<std::vector<std::string>>();
    absorbers_ = CreateAbsorbersFrom(absorber_names, rad);
  }

  // set rt solver
  if (my["rt-solver"]) {
    psolver = RTSolverFactory::CreateFrom(my["rt-solver"]);
  } else {
    psolver = nullptr;
  }

  char buf[80];
  snprintf(buf, sizeof(buf), "%.2f - %.2f", wmin_, wmax_);
  app->Log("Spectral range", buf);
  app->Log("Number of spectral bins", pgrid_->GetSize());
}

RadiationBand::~RadiationBand() {
  Application::Logger app("harp");
  app->Log("Destroy RadiationBand " + GetName());
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

void RadiationBand::WriteAsciiHeader(OutputParameters const *pout) const {
  if (!TestFlag(RadiationFlags::WriteBinRadiance)) return;

  char fname[80], number[6];
  snprintf(number, sizeof(number), "%05d", pout->file_number);
  snprintf(fname, sizeof(fname), "%s.radiance.%s.txt", GetName().c_str(),
           number);
  FILE *pfile = fopen(fname, "w");

  fprintf(pfile, "# Bin Radiances of Band %s: %.3g - %.3g\n", GetName().c_str(),
          wmin, wmax);
  fprintf(pfile, "# Ray output size: %lu\n", rayOutput_.size());

  fprintf(pfile, "# Polar angles: ");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%.3f", rad2deg(acos(rayOutput_[j].mu)));
  }
  fprintf(pfile, "\n");

  fprintf(pfile, "# Azimuthal angles: ");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%.3f", rad2deg(rayOutput_[j].phi));
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

  for (size_t i = 0; i < pgrid_->GetSize(); ++i) {
    fprintf(pfile, "%13.3g%12.3g", pgrid_->spec[i].wav1, pgrid_->spec[i].wav2);
    for (size_t j = 0; j < rayOutput_.size(); ++j) {
      fprintf(pfile, "%12.3f", toa_(i, j));
    }
    if (TestFlag(RadiationFlags::Normalize) && (wmax != wmin)) {
      fprintf(pfile, "%12.3g\n", pgrid_->spec[i].wght / (wmax - wmin));
    } else {
      fprintf(pfile, "%12.3g\n", pgrid_->spec[i].wght);
    }
  }

  fclose(pfile);
}

std::shared_ptr<RTSolver> RadiationBand::CreateRTSolveFrom(YAML::Node &my) {
  std::string rt_name_str = my["rt-solver"].as<std::string>();
  std::shared_ptr<RTSolver> psolver;

  if (rt_name_str == "Lambert") {
    psolver = std::make_shared<RTSolverLambert>(this);
#ifdef RT_DISORT
  } else if (rt_name_str == "Disort") {
    psolver = std::make_shared<RTSolverDisort>(this);
#endif
  } else {
    throw NotFoundError("RadiationBand", my["rt-solver"].as<std::string>());
  }

  return psolver;
}

RadiationBandContainer RadiationBandsFactory::CreateFrom(std::string filename) {
  Application::Logger app("harp");
  app->Log("Load Radiation bands from " + filename);

  std::vector<RadiationBandPtr> bands;

  std::ifstream stream(filename);
  if (stream.good() == false) {
    app->Error("Cannot open radiation bands file: " + filename);
  }
  YAML::Node node = YAML::Load(stream);

  if (!node["opacity-sources"]) {
    throw NotFoundError("LoadRadiationBand", "opacity-sources");
  }

  for (auto bname : node["bands"]) {
    auto p = std::make_shared<RadiationBand>(node, bname.as<std::string>());
    bands.push_back(p);
  }

  return bands;
}

RadiationBandContainer RadiationBandsFactory::CreateFrom(ParameterInput *pin,
                                                         std::string key) {
  std::vector<RadiationBandPtr> bands;

  auto rt_band_files =
      std::Vectorize<std::string>(pin->GetOrAddString("radiation", key, ""));

  for (auto &filename : rt_band_files) {
    auto &&tmp = CreateFrom(filename);
    bands.insert(bands.end(), tmp.begin(), tmp.end());
  }

  return bands;
}
