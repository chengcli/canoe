/** @file planets.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Mar 05, 2022 16:10:57 EST
 * @bug No known bugs.
 */

#ifndef PLANETS_HPP
#define PLANETS_HPP

//! planetocentric latitude to planetographic latitude
double centric2graphic(double clat_deg, double rerp);

//! planetographic latitude to planetocentric latitude
double graphic2centric(double glat_deg, double rerp);

//! planetographic latitude to planetocentric latitude for Jupiter
double jup_graphic2centric(double glat_deg);

//! planetocentric latitude to planetographic latitude for Jupiter
double jup_centric2graphic(double glat_deg);

#endif /* end of include guard PLANETS_HPP */
