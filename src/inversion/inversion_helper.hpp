
#include <string>
#include <vector>

// Eigen headers
#include <Eigen/Core>

class ParameterInput;
class MeshBlock;
class Inversion;

void read_observation_file(Eigen::VectorXd &target, Eigen::MatrixXd &icov,
                           std::string fname);

std::vector<Inversion *> create_inversion_queue(MeshBlock *pmb,
                                                ParameterInput *pin);

void delete_inversion_queue(std::vector<Inversion *> &fitq);

void gather_probability(std::vector<Inversion *> &fitq);
