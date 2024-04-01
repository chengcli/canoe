#ifndef SRC_PARAMETERIZATION_HPP_
#define SRC_PARAMETERIZATION_HPP_

class ParameterInput;
class MeshBlock;

namespace Parameterization {

void init_parameterization(MeshBlock *pmb, ParameterInput *pin);
void par_top_sponge_lyr_kl78(MeshBlock *pmb);

}  // namespace Parameterization

#endif  // SRC_PARAMETERIZATION_HPP_
