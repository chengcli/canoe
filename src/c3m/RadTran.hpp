// C/C++ headers
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// C3M headers
#include <configure.hpp>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Eigen::MatrixXd ReadVULCANPhotoAbsCrossSection(string Cross_Section_File);
Eigen::MatrixXd ReadStellarRadiationInput(string Solar_Input_File, double rad, double ref);

