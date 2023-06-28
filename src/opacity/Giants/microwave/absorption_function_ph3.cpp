/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to SNAP structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// ph3 data headers
#include "mwr_absorber_ph3.hpp"

// cli, 190801
inline double SQUARE(double x) { return x * x; }
inline double atm_from_bar(double P) { return 0.9869232667160128 * P; }
static double const pi = 3.14159265358979;
// end

double absorption_coefficient_PH3_radtran(double freq, double P, double T,
                                          double XH2, double XHe, double XPH3) {
  double dtaua = 0.0;

  double const PPH3 = atm_from_bar(P * XPH3);
  double const PH2 = atm_from_bar(P * XH2);
  double const PHe = atm_from_bar(P * XHe);

  // this flag controls whether to use the ad-hoc correction to the
  // PH3 absorption coefficients:   (TBD!)

  bool iph3cor = true;

  double const hc_over_kTref = 4.7959e-3;

  // sum the contributions of all lines:
  double const acon = 2.452e11;         // gives a in cm^-1
  double const fcon = 1.e-9 / 3.14159;  // old notes: 1.e12 ??
  double const texp = 3.5;

  for (int i = 0; i < NLINE_PH3; i++) {
    double frq0 = ph3_lines[i].FREQ * 1e-3;
    double intens = pow(10, ph3_lines[i].LGINT);
    double elow = ph3_lines[i].ELO;

    double tx = 300.0 / T;
    double fx = freq / frq0;
    double delnu = 1.419 * (PPH3 * 4.2154 + PH2 * 3.29283 + PHe * 1.6807);

    // don't use this form, for now:
    // double fline = fcon*fx*(delnu/(pow((freq-frq0),2)+ pow(delnu,2)) +
    //       delnu/(pow((freq+frq0),2) + pow(delnu,2)))

    double df = freq - frq0;
    double fline = fcon * fx * (delnu / (df * df + delnu * delnu));

    double tx35 = tx * tx * tx *
                  sqrt(tx);  // was pow(tx,texp), but gnu pow() is very slow
    double boltzf = exp(hc_over_kTref * elow * (1.0 - tx));
    dtaua += acon * PPH3 * intens * tx35 * boltzf * fx * fline;
  }

  if (iph3cor) dtaua *= 35.386 * pow(freq, -0.7969);

  return dtaua;  // in cm^-1
}

double absorption_coefficient_PH3_Hoffman(double freq, double P, double T,
                                          double XH2, double XHe, double XPH3) {
  double const factor =
      241.43;  // bar/P_cgs * 1/kb_cgs * (cm/nm)^2 / (300K) * (GHz/MHz)

  double dtaua = 0.0;

  double const PPH3 = P * XPH3;
  double const PH2 = P * XH2;
  double const PHe = P * XHe;

  double const hc_over_kTref = 4.7959e-3;

  // sum the contributions of all lines:
  for (int i = 0; i < NLINE_PH3; i++) {
    double frq0 = ph3_lines[i].FREQ * 1e-3;
    double intens = pow(10, ph3_lines[i].LGINT);
    double elow = ph3_lines[i].ELO;

    double tx = 300.0 / T;

    double gamma_H2, gamma_He, gamma_PH3;
    if (i < 40 && (std::abs(ph3_lines[i].QNp[1]) == 6 ||
                   (std::abs(ph3_lines[i].QNp[1]) == 3 &&
                    std::abs(ph3_lines[i].QNp[0]) < 8))) {
      intens *= 2.78;
      gamma_H2 = 1.4121;
      gamma_He = 0.7205;
      gamma_PH3 = 0.4976;
      // cout << "3x correction:  " << i << "  " << frq0 << endl;
    } else if (i < 40 &&
               (std::abs(ph3_lines[i].QNp[1]) == 3 &&
                std::abs(ph3_lines[i].QNp[0]) >= 8) &&
               std::abs(ph3_lines[i].QNp[0]) <= 26) {
      intens *= 36.65;
      gamma_H2 = 0.5978;
      gamma_He = 0.3050;
      gamma_PH3 = 3.1723;
      // cout << "37x correction:  " << i << "  " << frq0 << endl;
    } else {
      gamma_H2 = 3.2930;
      gamma_He = 1.6803;
      gamma_PH3 = 4.2157;
    }

    double fx = freq / frq0;
    double delnu = PPH3 * gamma_PH3 * tx +
                   pow(tx, 0.75) * (PH2 * gamma_H2 + PHe * gamma_He);

    double delnu_sq = SQUARE(delnu);
    double fline = fx / pi * delnu *
                   (1.0 / (SQUARE(freq - frq0) + delnu_sq) +
                    1.0 / (SQUARE(freq + frq0) + delnu_sq));

    double tx35 = tx * tx * tx * sqrt(tx);  // tx^3.5
    double boltzf = exp(hc_over_kTref * elow * (1.0 - tx));

    dtaua += factor * PPH3 * intens * tx35 * boltzf * fx * fline;
  }

  return dtaua;  // in cm^-1
}
