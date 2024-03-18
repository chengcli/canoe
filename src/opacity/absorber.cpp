// C/C++
#include <sstream>
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
  LoadOpacity(-1);
}

void Absorber::SetOpacityFile(std::string filename) {
  opacity_filename_ = filename;
}

void Absorber::LoadOpacity(int bid) {
  auto app = Application::GetInstance();
  auto log = app->GetMonitor("opacity");

  if (opacity_filename_.empty()) return;

  std::string full_path = app->FindResource(opacity_filename_);
  log->Log("Load opacity from " + full_path);
  LoadCoefficient(full_path, bid);
}

std::string Absorber::ToString() const {
  std::stringstream ss;
  ss << "Absorber: " << GetName();
  ss << "Opacity file: " << opacity_filename_;

  return ss.str();
}
