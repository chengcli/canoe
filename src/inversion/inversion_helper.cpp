// C/C++
#include <cstring>
#include <sstream>
#include <vector>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/vectorize.hpp>

// application
#include <application/application.hpp>

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
