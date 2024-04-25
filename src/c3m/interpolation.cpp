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

// Function to interpolate cross sections
Eigen::MatrixXd InterpolateCrossSection(Eigen::MatrixXd reference_wavelength, Eigen::MatrixXd input_wavelength, Eigen::MatrixXd input_cross_section){

int ref_size = reference_wavelength.size();
int input_size = input_wavelength.size();
VectorXd Output = VectorXd::Zero(ref_size);
double value;
int start = 0;
int end = input_size - 1;
int e1 = 0;
int pos = 0;

for(int ref_inx = 0; ref_inx < ref_size; ref_inx++){
  end = input_size - 1;
  start = 0;
//  Binary search
  while(end - start > 2){
    pos = (int)((start + end)/2);
// Element is in first half
    if((reference_wavelength(ref_inx) >= input_wavelength(start)) && (reference_wavelength(ref_inx) < input_wavelength(pos)))
    {
     end = pos;
    }
// Element is in second half
    if((reference_wavelength(ref_inx) <= input_wavelength(end)) && (reference_wavelength(ref_inx) >= input_wavelength(pos)))
    {start = pos;}
    if((reference_wavelength(ref_inx) < input_wavelength(start)))
    {end = start + 1;}
    if((reference_wavelength(ref_inx) > input_wavelength(end)))
    {start = end - 1;}}
// Interpolating the cross section
  value =  input_cross_section(start) + (reference_wavelength(ref_inx) - input_wavelength(start))*(input_cross_section(start) - input_cross_section(end))/(input_wavelength(start) - input_wavelength(end));
  
  if(value < 0){
  Output(ref_inx) = 0;}
  
  if(value >= 0){
  Output(ref_inx) = value;}
  

}

return Output;

}


Eigen::MatrixXd InterpolateQYield(Eigen::MatrixXd reference_wavelength, Eigen::MatrixXd input_wavelength, Eigen::MatrixXd qyield){
int ref_size = reference_wavelength.size();
int input_size = input_wavelength.size();
VectorXd Output = VectorXd::Zero(ref_size);
double value;
int start = 0;
int end = input_size - 1;
int e1 = 0;
int pos = 0;

for(int ref_inx = 0; ref_inx < ref_size; ref_inx++){
  end = input_size - 1;
  start = 0;
//  Binary search
  while(end - start > 2){
    pos = (int)((start + end)/2);
// Element is in first half
    if((reference_wavelength(ref_inx) >= input_wavelength(start)) && (reference_wavelength(ref_inx) < input_wavelength(pos)))
    {end = pos;}
  
// Element is in second half
    if((reference_wavelength(ref_inx) <= input_wavelength(end)) && (reference_wavelength(ref_inx) >= input_wavelength(pos)))
    {start = pos;}

    if((reference_wavelength(ref_inx) < input_wavelength(start)))
    {end = start + 1;}

    if((reference_wavelength(ref_inx) > input_wavelength(end)))
    {start = end - 1;}
  }
// Interpolating the qyield
  value =  qyield(start) + (reference_wavelength(ref_inx) - input_wavelength(start))*(qyield(start) - qyield(end))/(input_wavelength(start) - input_wavelength(end));
  
  if(value < 0){
  Output(ref_inx) = 0;}
  
  if((value >= 0) && (value < 1)){
  Output(ref_inx) = value;}
  
  if(value >= 1){
  Output(ref_inx) = 1;}
  

}

return Output;

}

