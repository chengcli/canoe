// C/C++ headers
#include <iostream>
#include <fstream>
// C3M headers
#include "PhotoChemistry.hpp"
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Athena++ header
#include <athena/parameter_input.hpp>

//NetCDF
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif


//Cantera header
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

//Print wavelengths
void printWavelength(Eigen::VectorXd wavelengths_)
{
  for (int i = 0; i < wavelengths_.size(); ++i) {
    std::cout << wavelengths_[i] << std::endl;
  }
}

//Print Cross Sections
void printCrossSection(Eigen::VectorXd crossSection_)
{

for (int i = 0; i < crossSection_.size(); ++i) {
    std::cout << crossSection_[i] << std::endl;
  }
}

// Calculate the photochemical reaction rate based on wavelength, cross section and actinic flux
double PhotoChemRate(Eigen::MatrixXd wavelengths_, Eigen::MatrixXd crossSection_, Eigen::MatrixXd Spectral_radiance)
{
int size = wavelengths_.size();
double sum = 0;
double sum1 = 0;
double pi = 3.141592;
double h = 6.626e-34;
double c = 3e8;
Spectral_radiance = (Spectral_radiance.array()*wavelengths_.array().transpose()/(h*c)).matrix();
MatrixXd integrd = (crossSection_.array()*Spectral_radiance.array()).matrix();
for(int  i = 0; i < size-1; i++){

  sum = sum + (integrd(i)*(wavelengths_(i+1) - wavelengths_(i)));
  sum1 = sum1 + (integrd(i+1)*(wavelengths_(i+1) - wavelengths_(i)));
}
sum = (sum + sum1)/2;
sum = sum*2*pi;
return sum;
}


double QPhotoChemRate(Eigen::MatrixXd wavelengths_, Eigen::MatrixXd d_wavelength, Eigen::MatrixXd crossSection_, Eigen::MatrixXd QYield, Eigen::MatrixXd Spectral_radiance)
{
int size = wavelengths_.size();
double sum = 0;
double sum1 = 0;
double pi = 3.141592;
double h = 6.626e-34;
double c = 3e8;
Spectral_radiance = (Spectral_radiance.array()*wavelengths_.array().transpose()/(h*c)).matrix();
MatrixXd integrd = (QYield.array()*crossSection_.array()*Spectral_radiance.array()).matrix();
sum = (integrd.array()*d_wavelength.array()).sum(); //Bhattacharya (5-1-23)

/*
for(int  i = 0; i < size-1; i++){

  sum = sum + (integrd(i)*(wavelengths_(i+1) - wavelengths_(i)));
  sum1 = sum1 + (integrd(i+1)*(wavelengths_(i+1) - wavelengths_(i)));
}
sum = (sum + sum1)/2;
*/

sum = sum*2*pi;
return sum;
}

// Function to read the photoionization cross sections from VULCAN database
// The output variable will contain wavelength (m) and photoionization cross section (m^2)
Eigen::MatrixXd  ReadVULCANPhotoIonCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID); 
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
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(photoion.c_str())*1E-4; //cm^2 to m^2
  
  num++;
  }
  InFile.close();
  
  return Output;
}


// Function to read the photodissociation cross sections from VULCAN database
// The output variable will contain wavelength (m) and photodissociation cross section (m^2)
Eigen::MatrixXd  ReadVULCANPhotoDissCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID); 
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
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(photodiss.c_str())*1E-4; //cm^2 to m^2
  
  num++;
  }
  InFile.close();
  return Output;
}

// Function to read cross section from VPL photochemical database
// The output variable contains the wavelength (m) and cross section (m^2)
Eigen::MatrixXd ReadAtmosCrossSection(string AtmosFileName){
  fstream InFile;
  InFile.open(AtmosFileName); 
  string wavlength;
  string photoabs;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  
  Eigen::MatrixXd Output(2, rows-4);
  InFile.open(AtmosFileName); 
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  while(InFile >> wavlength >> photoabs){
  
  Output(0, num) = atof(wavlength.c_str())*1E-10; //Angstrom to m
  Output(1, num) = atof(photoabs.c_str())*1E-4; //cm^2 to m^2

  
  num++;
  }
  InFile.close();
  return Output;
  
  
}

//Function to read cross section from Max Planck Mainz database
Eigen::MatrixXd ReadMPCrossSection(string MPFileName){
  fstream InFile;
  InFile.open(MPFileName);
  string wavlength;
  string photoabs;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();


  Eigen::MatrixXd Output(2, rows);
  InFile.open(MPFileName);
  while(InFile >> wavlength >> photoabs){

  Output(0, num) = atof(wavlength.c_str())*1E-10; //Angstrom to m 
  Output(1, num) = atof(photoabs.c_str())*1E-4; //cm^2 to m^2

  num++;
  }
  InFile.close();
  return Output;



}

//Function to read QY from VULCAN photochemistry database
Eigen::MatrixXd ReadQYield(string FileName){

fstream InFile;
  InFile.open(FileName); 
  //std::cout << FileName << std::endl;
  string wavlength, line;
  string photoabs;
  string photodiss;
  string photoion; 
  int num = 0;
  int rows = 0;
  int col = 0;
  while (getline(InFile, wavlength)){
  //std::cout << wavlength << std::endl;
  rows++;
  }
  InFile.close();
  InFile.open(FileName);
  getline(InFile, wavlength);
  getline(InFile, wavlength);
  while (getline(InFile, wavlength, ',')){
  col = col + 1;
  }
  
  InFile.close();
  //std::cout << "Col number " <<  (col-1)/(rows-2) << std::endl;
  InFile.open(FileName);
  Eigen::MatrixXd Output((col-1)/(rows-2)+1, rows-2);
  getline(InFile, wavlength);
  getline(InFile, wavlength);
  for(int ix =0; ix < rows-2; ix++){
  for(int ix2=0; ix2 <= (col-1)/(rows-2); ix2++){
  //std::cout << ix2 << " " << ix << std::endl; 
  if(ix2==0){
  getline(InFile, wavlength, ',');
  //std::cout << wavlength << std::endl;
  Output(0, ix) = atof(wavlength.c_str())*1E-9; //nm to m
  
  //std::cout << Output(0, ix) << std::endl;
  }
  if((ix2!=0) && (ix2 < (col-1)/(rows-2)) ){
  getline(InFile, wavlength, ',');
  //std::cout << wavlength << std::endl;
  Output(ix2, ix) = atof(wavlength.c_str());
  //std::cout << Output(ix2, ix) << std::endl;
  }
  if((ix2 == (col-1)/(rows-2))){
  getline(InFile, wavlength, '\n');
  //std::cout << wavlength << std::endl;
  Output(ix2, ix) = atof(wavlength.c_str());
  //std::cout << Output(ix2, ix) << std::endl;
  }
  
  
  
  }
  }
  InFile.close();
  /*
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(photodiss.c_str())*1E-4; //cm^2 to m^2
  
  num++;
  }
  InFile.close();
  */
  return Output;

}


//Function to read cross section from CalTech/JPL KINETIC7 corrected cross sections
//and converted to netCDF file
Eigen::MatrixXd ReadKINETICSCrossSection(int RxnIndex){
#if NETCDFOUTPUT
  int fileid, dimid, varid, err;
  string fname = "/home/ananyo/models/C3M/data/KINETICS/KINETICS7_Bhattacharya.nc";
  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  size_t nwaves;
  err = nc_inq_dimid(fileid, "wavelength", &dimid);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  err = nc_inq_dimlen(fileid, dimid, &nwaves);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  
  Eigen::VectorXd wavelength_(nwaves);
  Eigen::MatrixXd Output(2, nwaves);
  err = nc_inq_varid(fileid, "wavelength", &varid);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  err = nc_get_var_double(fileid, varid, wavelength_.data());
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  size_t nreactions = 571;

  Eigen::MatrixXd cross_section_(nreactions, nwaves);

  err = nc_inq_varid(fileid, "cross_section", &varid);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  err = nc_get_var_double(fileid, varid, cross_section_.data());

  for (size_t i = 0; i < wavelength_.size(); ++i){
    Output(0, i) = wavelength_(i)*1E-10; //Angstrom to m
    Output(1, i) =  cross_section_(RxnIndex,i)*1E-4; //cm^2 to m^2 
    //std::cout << Output(0, i) << " cross:  " << Output(1, i) << std::endl;
    }
  nc_close(fileid);

  return Output;
#endif

}



//Function to introduce custom opacity for a given atmosphere
Eigen::VectorXd handleCustomOpacity(string PlanetName,Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double Alt, Eigen::VectorXd wav){
int size = wav.size();
Eigen::VectorXd Output = VectorXd::Zero(size);

if(PlanetName == "Venus"){
   Output = VenusUV(NetworkName, Pres, Temp, Alt, wav);
  }


return Output;
}


//Function to add opacity due to unknown UV absorber on Venus
Eigen::VectorXd VenusUV(Cantera::ThermoPhase* NetworkName, double Pres, double Temp, double Alt, Eigen::VectorXd wav){


int size = wav.size();
Eigen::VectorXd COpacity = VectorXd::Zero(size);
//Altitude in km
if(Alt > 67){
  for(int i = 0; i < size; i++){
  COpacity(i) = 0.056*1E-3*exp( ((67E3-Alt)/3E3) - ((wav(i) - 3600E-10)/1000E-10));}
}

if((Alt <= 67) && (Alt > 58)){
  for(int i = 0; i < size; i++){
  COpacity(i) = 0.056*1E-3*exp( -1*((wav(i) - 3600E-10)/1000E-10));}
}


return COpacity;
}

// Function to read the Rayleigh scattering cross sections from VULCAN database
Eigen::MatrixXd ReadVULCANScatCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID);
  string wavlength;
  string scatcross;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();

  Eigen::MatrixXd Output(2, rows-1);
  InFile.open(VULCAN_ID);
  getline(InFile, wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,scatcross,'\n');


  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(scatcross.c_str())*1E-4; //cm^2 to m^2

  num++;
  }
  InFile.close();
  return Output;
}

