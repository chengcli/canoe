// Athena++ headers
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>

// canoe headers
#include "face_reconstruct.hpp"

FaceReconstruct::FaceReconstruct(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {}
