/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to SNAP structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// cli, 190801
inline double SQUARE(double x) { return x * x; }
// end

double absorption_coefficient_cloud(double freq, double P, double T,
                                    double rho_NH3_H2O, double rho_H2O,
                                    double rho_NH4SH, double rho_NH3,
                                    double cfliq, double cfwice,
                                    double cfaice) {
  //  calculate absorption by clouds
  //
  //  Unless otherwise stated, equations are from:
  //  Atmospheric Remote Sensing by Microwave Radiometer, Ed. Michael A.
  //  Janssen, Wiley, Page 332. 1993.

  // liquid cloud (assumed pure H2O, ref: Grody, N.C.):
  // (valid for FREQ < 100 GHz)
  double nu0 = 160.0 * exp(7.21 * (1.0 - 287. / T));
  double acldl = 1.24e-1 * freq * freq * nu0 / (freq * freq + nu0 * nu0);

  // H2O ice cloud (ref: Gasiewski, A.J.):
  double tt = 300.0 / T - 1.0;
  double tt1 = 1.0 + tt;
  double tt2 = 0.0073 + tt;
  double aa = (0.00504 + 0.0062 * tt) * exp(-22.1 * tt);
  double bb =
      1.0e-4 * (0.502 - 0.131 * tt) / tt1 + 0.542e-6 * tt1 * tt1 / (tt2 * tt2);
  double ei = aa / freq + bb * freq;
  double er = 3.1884 + 9.1e-4 * (T - 273.0) + 2.0;

  double acldw = 1.885e0 * freq * ei / (er * er + ei * ei);

  // NH3 ice cloud (approximation by S.Gulkis, private notes, 2007):
  // (for simplicity, NH4SH is lumped together with this, for now)
  double aclda = 1.13e-4 * SQUARE(freq) / (42.25 + SQUARE(6e-5 * freq));
  // old expression:
  // double aclda = 0.27 * freq * freq;

  // total cloud absorption:
  return cfliq * acldl * rho_NH3_H2O + cfwice * acldw * rho_H2O +
         cfaice * aclda * (rho_NH4SH + rho_NH3);
}
