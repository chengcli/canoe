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
#include <cantera/kinetics/Reaction.h>
#include <cantera/kinetics/ReactionData.h>
#include <cantera/kinetics/ReactionRate.h>
#include <cantera/kinetics/MultiRate.h>
#include <cantera/base/Solution.h>
#include <cantera/thermo.h>

// C3M headers
#include <configure.hpp>
#include "CustomRate.hpp"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


void handleCustomChemistry(string PlanetName, Cantera::Kinetics* NetworkName, double Pres, double Temp){

  if(PlanetName == "JupiterAurora")
    {JupiterAuroralChemistry_VT(NetworkName, Pres, Temp);
    }


}


void JupiterAuroralChemistry_VT(Cantera::Kinetics* NetworkName, double Pres, double Temp){
  int numReaction = NetworkName->nReactions();
for(int inumRxn =  0; inumRxn < numReaction; inumRxn++){
  auto& react = *(NetworkName->reaction(inumRxn));
  std::string Equation = react.equation();
  //int pos = Equation.find("=");
  //Equation.replace(pos, 2, "->");
 /*

2 H + H2 => H2_v1 + H2
2 H + H2 => H2_v2 + H2
2 H + H2 => H2_v3 + H2
2 H + H2 => H2_v4 + H2
2 H + H2 => H2_v5 + H2
2 H + H2 => H2_v6 + H2
2 H + H2 => H2_v7 + H2
2 H + H2 => H2_v8 + H2
H+ + H2_v4 => H + H2+
H+ + H2_v5 => H + H2+
H+ + H2_v6 => H + H2+
H+ + H2_v7 => H + H2+
H+ + H2_v8 => H + H2+
E + H3+ => H + H2_v1
E + H3+ => H + H2_v2
E + H3+ => H + H2_v3
E + H3+ => H + H2_v4
E + H3+ => H + H2_v5
E + H3+ => H + H2_v6
E + H3+ => H + H2_v7
E + H3+ => H + H2_v8
H + H2_v1 => H + H2
H + H2_v2 => H + H2_v1
H + H2_v3 => H + H2_v2
H + H2_v4 => H + H2_v3
H + H2_v5 => H + H2_v4
H + H2_v6 => H + H2_v5
H + H2_v7 => H + H2_v6
H + H2_v8 => H + H2_v7
H+ + H2_v1 => H+ + H2
H+ + H2_v2 => H+ + H2_v1
H+ + H2_v3 => H+ + H2_v2
H+ + H2_v4 => H+ + H2_v3
H+ + H2_v5 => H+ + H2_v4
H+ + H2_v6 => H+ + H2_v5
H+ + H2_v7 => H+ + H2_v6
H+ + H2_v8 => H+ + H2_v7
H2 + H2_v1 => 2 H2
H2 + H2_v2 => H2 + H2_v1
H2 + H2_v3 => H2 + H2_v2
H2 + H2_v4 => H2 + H2_v3
H2 + H2_v5 => H2 + H2_v4
H2 + H2_v6 => H2 + H2_v5
H2 + H2_v7 => H2 + H2_v6
H2 + H2_v8 => H2 + H2_v7
H2 + H2_v1 => H2 + H2_v1
2 H2_v1 => H2 + H2_v2
H2_v1 + H2_v2 => H2 + H2_v3
H2_v1 + H2_v3 => H2 + H2_v4
H2_v1 + H2_v4 => H2 + H2_v5
H2_v1 + H2_v5 => H2 + H2_v6


*/


//Add the reactions and corresponding reaction rate expressions
//VT chemistry for H2-H2 collisions based on Cravens et al., 1987

  if(Equation == "H2 + H2_v1 => 2 H2"){
   int vVal = 1;
   double rateValue = 1.0;
   //std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v2 => H2 + H2_v1"){
   int vVal = 2;
   double rateValue = 1.0;
  // std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v3 => H2 + H2_v2"){
   int vVal = 3;
   double rateValue = 1.0;
  // std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v4 => H2 + H2_v3"){
   int vVal = 4;
   double rateValue = 1.0;
 //  std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v5 => H2 + H2_v4"){
   int vVal = 5;
   double rateValue = 1.0;
  // std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v6 => H2 + H2_v5"){
   int vVal = 6;
   double rateValue = 1.0;
  // std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v7 => H2 + H2_v6"){
   int vVal = 7;
   double rateValue = 1.0;
 //  std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H2 + H2_v8 => H2 + H2_v7"){
   int vVal = 8;
   double rateValue = 1.0;
 //  std::cout << "calculation is happening" << std::endl;
   if(Temp < 500){
   rateValue = 1.59E-11*exp(-78.75/pow(Temp, 1/3));
   }
   if(Temp >= 500){
   rateValue = 2.14E-9*exp(-117.7/pow(Temp, 1/3));
   }

   double del = 1.3 - (0.3*(Temp - 500)/700);
   rateValue = rateValue*exp(del*(vVal -1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }
//VT chemistry for H-H2 collisions based on Cravens et al., 1987

  if(Equation == "H + H2_v1 => H + H2"){
   int vVal = 1;
 //  std::cout << "calculation is happening" << std::endl;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H + H2_v2 => H + H2_v1"){
   int vVal = 2;
 //  std::cout << "calculation is happening" << std::endl;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }
 
  if(Equation == "H + H2_v3 => H + H2_v2"){
   int vVal = 3;
 //  std::cout << "calculation is happening" << std::endl;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H + H2_v4 => H + H2_v3"){
   int vVal = 4;
//   std::cout << "calculation is happening" << std::endl;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H + H2_v5 => H + H2_v4"){
   int vVal = 5;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H + H2_v6 => H + H2_v5"){
   int vVal = 6;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

  if(Equation == "H + H2_v7 => H + H2_v6"){
   int vVal = 7;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "H + H2_v8 => H + H2_v7"){
   int vVal = 8;
   double rateValue = 1.4E-13*exp( (Temp/125) - pow((Temp/577), 2));
   rateValue = rateValue*exp(0.714*(vVal - 1));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }
//VV chemistry for H2-H2 collisions based on Cravens et al., 1987 
 if(Equation == "H2_v1 + H2 => H2 + H2_v1"){
   double rateValue = 9.12E-12*exp(-24/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "2 H2_v1 => H2 + H2_v2"){
   double rateValue = 6.83E-12*exp(-15.58/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "H2_v1 + H2_v2 => H2 + H2_v3"){
   double rateValue = 2.46E-11*exp(-21.84/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "H2_v1 + H2_v3 => H2 + H2_v4"){
   double rateValue = 1.65E-11*exp(-20.58/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "H2_v1 + H2_v4 => H2 + H2_v5"){
   double rateValue = 6.43E-11*exp(-27.54/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  }

 if(Equation == "H2_v1 + H2_v5 => H2 + H2_v6"){
   double rateValue = 1.01E-10*exp(-32.69/pow(Temp, 1/3));
   NetworkName->setMultiplier(inumRxn, rateValue);
  } 
 
 

}}
