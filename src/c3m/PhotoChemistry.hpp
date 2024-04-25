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

// C3M headers
#include <configure.hpp>
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


void printWavelength(Eigen::VectorXd wavelengths_);
void printCrossSection(Eigen::VectorXd crossSection_);
double PhotoChemRate(Eigen::MatrixXd wavelengths_, Eigen::MatrixXd crossSection_, Eigen::MatrixXd Spectral_radiance);
double QPhotoChemRate(Eigen::MatrixXd wavelengths_, Eigen::MatrixXd d_wavelength, Eigen::MatrixXd crossSection_, Eigen::MatrixXd QYield, Eigen::MatrixXd Spectral_radiance);
Eigen::MatrixXd ReadJPLCrossSection(string JPL_FileName, string Reactant);
Eigen::MatrixXd ReadVULCANPhotoIonCrossSection(string VULCAN_ID);
Eigen::MatrixXd ReadVULCANPhotoDissCrossSection(string VULCAN_ID);
Eigen::MatrixXd ReadVULCANScatCrossSection(string VULCAN_ID);
Eigen::MatrixXd ReadAtmosCrossSection(string AtmosFileName);
Eigen::MatrixXd ReadMPCrossSection(string MPFileName);
Eigen::MatrixXd ReadKINETICSCrossSection(int RxnIndex);
Eigen::MatrixXd ReadQYield(string FileName);


Eigen::VectorXd handleCustomOpacity(string PlanetName, Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double Alt, Eigen::VectorXd wav);
Eigen::VectorXd VenusUV(Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double Alt, Eigen::VectorXd wav);

