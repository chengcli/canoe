/** @file planets.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Mar 05, 2022 16:13:16 EST
 * @bug No known bugs.
 */

// climath headers
#include <climath/core.h>

// harp2 headers
#include "planets.hpp"

double centric2graphic(double clat_deg, double rerp) {
  if (std::abs(clat_deg) == 90.)
    return clat_deg;
  else
    return rad2deg(atan(tan(deg2rad(clat_deg))*(rerp*rerp)));
}

double graphic2centric(double glat_deg, double rerp) {
  if (std::abs(glat_deg) == 90.)
    return glat_deg;
  else
    return rad2deg(atan(tan(deg2rad(glat_deg))/(rerp*rerp)));
}

double jup_graphic2centric(double glat_deg) {
  double jup_rerp = 1.07;
  return graphic2centric(glat_deg, jup_rerp);
}

double jup_centric2graphic(double clat_deg) {
  double jup_rerp = 1.07;
  return centric2graphic(clat_deg, jup_rerp);
}
