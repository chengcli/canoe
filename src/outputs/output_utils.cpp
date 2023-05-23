#include "output_utils.hpp"

#include "user_outputs.hpp"

extern DiagnosticTable diag_table;

int get_num_variables(std::string grid, AthenaArray<Real> const& data) {
  int nvar;
  if (grid == "--C" || grid == "--F") {
    nvar = data.GetDim2();
  } else if (grid == "---") {
    nvar = data.GetDim1();
  } else {
    nvar = data.GetDim4();
  }

  return nvar;
}

std::string get_grid_type(std::string name) {
  int nouts = diag_table.size();
  for (int i = 0; i < nouts; ++i) {
    if (diag_table[i][0] == name) {
      return diag_table[i][3];
    }
  }

  return "";
}

std::string get_units(std::string name) {
  int nouts = diag_table.size();
  for (int i = 0; i < nouts; ++i) {
    if (diag_table[i][0] == name) {
      return diag_table[i][2];
    }
  }

  return "";
}

std::string get_long_name(std::string name) {
  int nouts = diag_table.size();
  for (int i = 0; i < nouts; ++i) {
    if (diag_table[i][0] == name) {
      return diag_table[i][1];
    }
  }

  return "";
}
