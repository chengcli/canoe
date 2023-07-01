#ifndef SRC_SETUP_MESH_HPP_
#define SRC_SETUP_MESH_HPP_

class ParameterInput;
class Mesh;

void setup_mesh(ParameterInput*& pin, Mesh*& mesh);
void destroy_mesh(ParameterInput*& pin, Mesh*& mesh);

#endif
