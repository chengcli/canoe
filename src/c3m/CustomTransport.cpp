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
#include "CustomTransport.hpp"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


VectorXd handleCustomMolecularDiffusion(string PlanetName, Cantera::ThermoPhase* NetworkName, double Pres, double Temp, VectorXd mWt){
 int nsp = NetworkName->nSpecies();
 VectorXd DiffConst = VectorXd::Zero(nsp);

  if(PlanetName == "JupiterAurora")
    { 
      DiffConst = JupiterMolDiff(NetworkName, Pres, Temp, mWt);
    }
 return DiffConst;
}


VectorXd JupiterMolDiff(Cantera::ThermoPhase* NetworkName, double Pres, double Temp, VectorXd mWt){
 int nsp = NetworkName->nSpecies();
 VectorXd DiffConst = VectorXd::Zero(nsp);
 VectorXd col_freq = VectorXd::Zero(nsp);
 VectorXd mol_fr = VectorXd::Zero(nsp);
 VectorXd I = VectorXd::Ones(nsp);
 int elementIndex = 0; //NetworkName->elementIndex("H2");
 NetworkName->getMoleFractions(&mol_fr[0]);
 double r = 2.7E-10; //Collision radius of H2 (m)
//Get total number density of H2
 double totalDensity = NetworkName->molarDensity(); //kmol/m^3
 double nH2 = totalDensity*1E3*6.022E23*mol_fr(elementIndex); //kmol/m^3 -> #/m^3
//Compute collision frequency
 col_freq = 6.022e23*2*3.14*1.38E-23*Temp*(mWt + (I*mWt(elementIndex)))/(mWt(elementIndex)*1E-3);
 col_freq = (col_freq.array()/mWt.array()).sqrt().matrix();
//Compute binary diffusion coefficient
 col_freq = col_freq*2*nH2*r*r;
 DiffConst = (1/col_freq.array()).matrix();
 DiffConst = (DiffConst.array()/mWt.array()).matrix()*1E3*1.38E-23*6.022e23*Temp;
// DiffConst(1) = 3.13E-5*7.339E11*pow(Temp, 0.765);
 DiffConst(0) = 0.0;
// std::cout << "Diffusion constant" << std::endl;
// std::cout << DiffConst << std::endl;
 return DiffConst;
}

