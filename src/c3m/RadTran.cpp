// C/C++ headers
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// C3M headers
#include <configure.hpp>

#include "Vapors.hpp"
using Eigen::MatrixXd;
using Eigen::VectorXd;


// Function to read photoabsorption cross section from VULCAN database
// The output matrix will contain wavelength (nm) and photoabsorption cross section (cm^2)
Eigen::MatrixXd  ReadVULCANPhotoAbsCrossSection(string Cross_Section_File){
  fstream InFile;
  InFile.open(Cross_Section_File); 
  string wavlength;
  string photoabs;
  string photodiss;
  string photoion; 
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  Eigen::MatrixXd Output(2, rows-1);
  InFile.open(Cross_Section_File); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm -> m
  Output(1, num) = atof(photoabs.c_str())*1E-4; // cm^2 -> m^2
  
  num++;
  }
  
  InFile.close();
  return Output;
}

// Function to read stellar radiation input
// The output matrix will contain reference wavelength (unit) and spectral irradiance (unit)
Eigen::MatrixXd ReadStellarRadiationInput(string Solar_Input_File, double rad, double ref){
  fstream InFile;
  InFile.open(Solar_Input_File); 
  string wavlength;
  string irradiance;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  Eigen::MatrixXd Output(2, rows-1);
  InFile.open(Solar_Input_File); 
  getline(InFile,wavlength);
  while(InFile >> wavlength >> irradiance){
    
  Output(0, num) = stof(wavlength.c_str())*1E-10; //angstrom -> m
  Output(1, num) = stof(irradiance.c_str())*1E10*ref*ref/(rad*rad); //Conversion to SI units (W/m^3)
  //std::cout << irradiance << " " << Output(1, num) << std::endl;
  num++;
  }
  InFile.close();
  return Output;


}

