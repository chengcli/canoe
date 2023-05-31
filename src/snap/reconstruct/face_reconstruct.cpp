// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// snap
#include "face_reconstruct.hpp"

FaceReconstruct::FaceReconstruct(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {}
