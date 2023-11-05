// C/C++
#include <string>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// opacity
#include "absorber.hpp"

Absorber::Absorber(std::string name) : NamedGroup(name), opacity_filename_("") {
  Application::Logger app("opacity");
  app->Log("Create Absorber " + name);
}

Absorber::~Absorber() {
  Application::Logger app("opacity");
  app->Log("Destroy Absorber " + GetName());
}

void Absorber::LoadOpacityFromFile(std::string filename) {
  SetOpacityFile(filename);
  LoadOpacity();
}

void Absorber::SetOpacityFile(std::string filename) {
  opacity_filename_ = filename;
}

void Absorber::LoadOpacity() {
  auto app = Application::GetInstance();
  auto log = app->GetMonitor("opacity");

  if (opacity_filename_.empty()) return;

  try {
    std::string full_path = app->FindInputFile(opacity_filename_);
    log->Log("Load opacity from " + full_path);
    LoadCoefficient(full_path, 0);
  } catch (NotFoundError const& e) {
    std::stringstream ss;
    ss << e.what() << std::endl;
    log->Warn(ss.str());
  }
}
