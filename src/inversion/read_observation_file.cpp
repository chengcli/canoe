// C/C++ headers
#include <athena/athena.hpp>
#include <cstring>
#include <sstream>
#include <utils/fileio.hpp>

#include "inversion.hpp"

#define MAX_LINE 512
void read_observation_file(Eigen::VectorXd &target, Eigen::MatrixXd &icov,
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
  target.resize(rows);
  icov.resize(rows, rows);

  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    sscanf(pl, "%lf", &target(i));
  }

  // inverse covariance matrix
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok_r(pl, " ");
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &icov(i, j));
      p = strtok_r(NULL, " ");
    }
  }

  fclose(fp);
}
#undef MAX_LINE
