// C/C++ headers
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

// Eigen
#include <Eigen/Core>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"
#include "profile_inversion.hpp"

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  Application::Logger app("inversion");
  app->Log("Initialize Inversion");
}

std::vector<InversionPtr> InversionsFactory::CreateFrom(MeshBlock *pmb,
                                                        ParameterInput *pin) {
  Application::Logger app("inversion");
  app->Log("Create inversion queue");

  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names =
      Vectorize<std::string>(str.c_str(), " ,");

  std::vector<InversionPtr> all_fits;
  InversionPtr pfit;

  int jl = pmb->js, ju = pmb->js;
  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = std::make_shared<VLAProfileInversion>(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = std::make_shared<JunoProfileInversion>(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      throw NotFoundError("CreateAllInversions", p);
    }

    ju += pfit->GetSteps();
    pfit->SetStepRange(jl, ju);
    jl = ju;
    all_fits.push_back(std::move(pfit));
  }

  app->Log("Number of inversions = " + std::to_string(all_fits.size()));

  return all_fits;
}

namespace InversionHelper {
#define MAX_LINE 512
void read_observation_file(Eigen::VectorXd *target, Eigen::MatrixXd *icov,
                           std::string fname) {
  std::stringstream msg;
  FILE *fp = fopen(fname.c_str(), "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in ProfileInversion::ReadObseravtionFile"
        << std::endl
        << fname << " cannot be opened.";
    ATHENA_ERROR(msg);
  }
  char line[MAX_LINE], *pl;

  int rows;
  // header
  pl = NextLine(line, MAX_LINE, fp);

  // target values
  sscanf(pl, "%d", &rows);
  target->resize(rows);
  icov->resize(rows, rows);

  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    sscanf(pl, "%lf", &(*target)(i));
  }

  // inverse covariance matrix
  char *saveptr;
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok_r(pl, " ", &saveptr);
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &(*icov)(i, j));
      p = strtok_r(NULL, " ", &saveptr);
    }
  }

  fclose(fp);
}
#undef MAX_LINE

void gather_probability(std::vector<Inversion *> const &fitq) {
  auto qlast = fitq[fitq.size() - 1];

  // replace the log probability by the last one
  for (auto q : fitq) {
    int nwalker = q->GetWalkers();

    for (int k = 0; k < nwalker; ++k) {
      q->SetLogProbability(k, qlast->GetLogProbability(k));
    }
  }
}
}  // namespace InversionHelper
