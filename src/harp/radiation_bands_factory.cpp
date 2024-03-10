// external
#include <application/application.hpp>

// utils
#include <utils/vectorize.hpp>

// harp
#include "radiation_band.hpp"

std::map<std::string, int> RadiationBandsFactory::band_id_ = {};
int RadiationBandsFactory::last_band_id_ = 0;

RadiationBandContainer RadiationBandsFactory::CreateFrom(std::string filename) {
  Application::Logger app("harp");
  app->Log("Load Radiation bands from " + filename);

  std::vector<RadiationBandPtr> bands;

  std::ifstream stream(filename);
  if (stream.good() == false) {
    app->Error("Cannot open radiation bands file: " + filename);
  }
  YAML::Node rad = YAML::Load(stream);

  for (auto bname : rad["bands"]) {
    band_id_[bname.as<std::string>()] = last_band_id_++;

    auto p = std::make_shared<RadiationBand>(bname.as<std::string>(), rad);
    bands.push_back(p);
  }

  return bands;
}

RadiationBandContainer RadiationBandsFactory::CreateFrom(ParameterInput *pin,
                                                         std::string key) {
  std::vector<RadiationBandPtr> bands;

  auto rt_band_files =
      Vectorize<std::string>(pin->GetOrAddString("radiation", key, "").c_str());

  for (auto &filename : rt_band_files) {
    auto &&tmp = CreateFrom(filename);
    bands.insert(bands.end(), tmp.begin(), tmp.end());
  }

  return bands;
}
