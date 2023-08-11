// C/C++
#include <string>

// application
#include <application/application.hpp>

// harp
#include <harp/radiation_band.hpp>

// utils
#include <utils/extract_substring.hpp>
#include <utils/fileio.hpp>
#include <utils/parameter_map.hpp>

// opacity
#include <opacity/Giants/microwave/mwr_absorbers.hpp>
#include <opacity/Giants/freedman.hpp>
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/correlatedk_absorber.hpp>
#include <opacity/hitran_absorber.hpp>


void RadiationBand::addAbsorberInfrared(YAML::Node &node) {
  auto name = node["name"].as<std::string>();
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  std::string data_file;

  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (node["data"]) {
    data_file = node["data"].as<std::string>();
  } else {
    data_file = "kcoeff." + myfile_ + "-" + name_ + ".nc";
  }

  auto app = Application::GetInstance();
  std::string full_path;
  try {
    full_path = app->FindInputFile(data_file);
  } catch (NotFoundError const &e) {
    auto log = app->GetMonitor("harp");
    std::stringstream ss;
    ss << e.what() << std::endl;
    ss << name << " will not be loaded." << std::endl;
    log->Warn(ss.str());
  }

  if (name == "H2-H2-CIA") {
    auto ab = std::make_unique<XizH2H2CIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2-He-CIA") {
    auto ab = std::make_unique<XizH2HeCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CH4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H2") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H6") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "freedman_simple") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);
    absorbers_.push_back(std::move(ab));

  } else if (name == "freedman_mean") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);
    absorbers_.push_back(std::move(ab));

  } else if (strncmp(name.c_str(), "ck-", 3) == 0) {
    auto aname = extract_second(name, "-");
    auto ab = std::make_unique<CorrelatedKAbsorber>(aname, species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberInfrared", "Absorber " + name);
  }
}

// radiation
void RadiationBand::addAbsorberRadio(YAML::Node& node) {
  auto name = node["name"].as<std::string>();
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "NH3") {
    auto ab = std::make_unique<gp::MwrAbsorberNH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2O") {
    auto ab = std::make_unique<gp::MwrAbsorberH2O>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2S") {
    auto ab = std::make_unique<gp::MwrAbsorberH2S>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "PH3") {
    auto ab = std::make_unique<gp::MwrAbsorberPH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CIA") {
    auto ab = std::make_unique<gp::MwrAbsorberCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "Electron") {
    auto ab = std::make_unique<gp::MwrAbsorberElectron>(species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberRadio", "Absorber " + name);
  }
}
