// C/C++
#include <memory>
#include <string>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"

CompositionInversion::CompositionInversion(MeshBlock *pmb,
                                           YAML::Node const &node)
    : Inversion(pmb, "composition") {
  Application::Logger app("inversion");
  app->Log("Initializing CompositionInversion");
  char buf[80];

  // species id
  idx_ =
      Vectorize<int>(pin->GetString("inversion", name + ".variables").c_str());

  // read in prior
  for (auto m : idx_) {
    if (m == IDN) {  // change temperature
      Xstd_[IDN] = pin->GetReal("inversion", name + ".tem.std");
      app->Log(name + "::temperature std = " + std::to_string(Xstd_[IDN]));
    } else {
      Xstd_[m] = pin->GetReal("inversion", name + ".qvapor" +
                                               std::to_string(m) + ".std.gkg") /
                 1.E3;
      snprintf(buf, sizeof(buf), "%s::vapor %d standard deviation",
               name.c_str(), m);
      app->Log(buf + std::to_string(Xstd_[m]));
    }
  }
}

void CompositionInversion::UpdateModel(std::vector<Real> const &par) const {
  Application::Logger app("inversion");
  app->Log("UpdateConcentration");

  // int is = pblock_->is, ie = pblock_->ie;
}
