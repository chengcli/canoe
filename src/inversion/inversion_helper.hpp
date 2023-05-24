#ifndef SRC_INVERSION_INVERSION_HELPER_HPP_ 
#define SRC_INVERSION_INVERSION_HELPER_HPP_

#include <string>
#include <vector>

// Eigen headers
#include <Eigen/Core>

class ParameterInput;
class MeshBlock;
class Inversion;

void read_observation_file(Eigen::VectorXd *target, Eigen::MatrixXd *icov,
                           std::string fname);

std::vector<Inversion *> create_inversion_queue(MeshBlock *pmb,
                                                ParameterInput *pin);

void gather_probability(std::vector<Inversion *> const& fitq);

#endif  // SRC_INVERSION_INVERSION_HELPER_HPP_
