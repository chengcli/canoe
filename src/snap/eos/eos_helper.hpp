#include <athena/defs.hpp>

class MeshBlock;
template <typename T>
class AthenaArray;

void apply_vapor_limiter(MeshBlock *pmb, AthenaArray<Real> &u);
void check_hydro_variables(MeshBlock *pmb, AthenaArray<Real> &u);
