#ifndef SRC_SNAP_EOS_EOS_HELPER_HPP_
#define SRC_SNAP_EOS_EOS_HELPER_HPP_

#include <athena/defs.hpp>

class MeshBlock;
template <typename T>
class AthenaArray;

void apply_vapor_limiter(AthenaArray<Real> *pu, MeshBlock *pmb);
void check_hydro_variables(AthenaArray<Real> *pu, MeshBlock *pmb);

#endif  // SRC_SNAP_EOS_EOS_HELPER_HPP_
