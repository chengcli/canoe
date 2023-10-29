//! \file dirty.hpp
//! \brief dirty.hpp contains dirty functions that should go somewhere else

#ifndef SRC_DIRTY_HPP_
#define SRC_DIRTY_HPP_

// helper functions, will be moved in the future
int find_pressure_level_lesser(Real pres, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie) {
  for (int i = is; i <= ie; ++i)
    if (w(IPR, k, j, i) < pres) return i;

  return ie + 1;
}

#endif  // SRC_DIRTY_HPP_
