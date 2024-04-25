// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// C/C++ headers
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// Athena++ header
#include <athena/parameter_input.hpp>

// Cantera headers
#include <cantera/kinetics/Kinetics.h>
#include <cantera/base/ct_defs.h>
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>
#include <cantera/kinetics/MultiRate.h>
#include <cantera/base/Solution.h>
#include <cantera/thermo.h>

// C3M headers
#include <configure.hpp>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd handleCustomMolecularDiffusion(string PlanetName, Cantera::ThermoPhase* NetworkName, double Pres, double Temp, VectorXd mWt);

VectorXd JupiterMolDiff(Cantera::ThermoPhase* NetworkName, double Pres, double Temp, VectorXd mWt);
