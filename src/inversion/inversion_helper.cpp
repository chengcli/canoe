// C/C++
#include <cstring>
#include <sstream>
#include <vector>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"
#include "inversion_helper.hpp"
#include "profile_inversion.hpp"

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
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok_r(pl, " ");
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &(*icov)(i, j));
      p = strtok_r(NULL, " ");
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

std::vector<Inversion *> create_inversion_queue(MeshBlock *pmb,
                                                ParameterInput *pin) {
  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names =
      Vectorize<std::string>(str.c_str(), " ,");

  std::vector<Inversion *> fitq;

  Inversion *pfit;
  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = new VLAProfileInversion(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = new JunoProfileInversion(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      Debugger::Fatal("new_inversion_queue", "task::" + p, "unrecognized");
    }
    fitq.push_back(pfit);
  }

  int jl = pmb->js;
  for (auto q : fitq) {
    q->InitializePositions();
    q->setX2Indices(jl);
    jl += q->getX2Span();
  }

  return pfit;
}
