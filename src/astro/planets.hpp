#ifndef SRC_ASTRO_PLANETS_HPP_
#define SRC_ASTRO_PLANETS_HPP_

//! planetocentric latitude to planetographic latitude
double centric2graphic(double clat_deg, double rerp);

//! planetographic latitude to planetocentric latitude
double graphic2centric(double glat_deg, double rerp);

//! planetographic latitude to planetocentric latitude for Jupiter
double jup_graphic2centric(double glat_deg);

//! planetocentric latitude to planetographic latitude for Jupiter
double jup_centric2graphic(double glat_deg);

#endif  // SRC_ASTRO_PLANETS_HPP_
