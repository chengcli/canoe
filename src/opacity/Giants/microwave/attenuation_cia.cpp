/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to canoe structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// CIA data headers
#include "mwr_absorber_cia.hpp"

// cli, 190801
inline double SQUARE(double x) { return x * x; }
inline double Pa_from_bar(double P) { return 1.0e5 * P; }
static double const kBoltz_mks = 1.3806504e-23;  // J/K
static double const Lo_mks =
    2.68719e25;  // Loschmidt number molecules/m^3 at STP
static double const c_cgs = 2.99792458e+10;  // speed of light
// end

double absorption_coefficient_CIA(double freq, double P, double T, double XH2,
                                  double XHe, double XCH4, double mix) {
  double const factor = Pa_from_bar(P) / (kBoltz_mks * T) / Lo_mks;
  double const amagat_he = XHe * factor;
  double const amagat_h2 = XH2 * factor;
  double const amagat_ch4 = XCH4 * factor;

  double nu_invcm = freq * 1e9 / c_cgs;
  double logT = log(T);
  int Inu0 = CIA_Orton::find_index_wavenumber(nu_invcm);
  int It0 = CIA_Orton::find_index_temp(logT);

  // if we are below Orton's spectral range, assume a nu^2 dependence
  double nu_scale = 1.0;
  if (Inu0 < 0) {
    Inu0 = 0;
    nu_scale = SQUARE(nu_invcm / CIA_Orton::wn[0]);
  }

  // mix = 0 (equilibrium), 1 (normal, high temp)
  double h2h2, h2he, h2ch4;
  if (mix == 0.0) {
    h2h2 = amagat_h2 * amagat_h2 *
           exp(CIA_Orton::CIA_H2_H2_e(It0, Inu0, nu_invcm, logT));
    h2he = amagat_h2 * amagat_he *
           exp(CIA_Orton::CIA_H2_He_e(It0, Inu0, nu_invcm, logT));
    h2ch4 = amagat_h2 * amagat_ch4 *
            exp(CIA_Orton::CIA_H2_CH4_e(It0, Inu0, nu_invcm, logT));
  } else if (mix == 1.0) {
    h2h2 = amagat_h2 * amagat_h2 *
           exp(CIA_Orton::CIA_H2_H2_n(It0, Inu0, nu_invcm, logT));
    h2he = amagat_h2 * amagat_he *
           exp(CIA_Orton::CIA_H2_He_n(It0, Inu0, nu_invcm, logT));
    h2ch4 = amagat_h2 * amagat_ch4 *
            exp(CIA_Orton::CIA_H2_CH4_n(It0, Inu0, nu_invcm, logT));
  } else {
    h2h2 = amagat_h2 * amagat_h2 *
           (mix * exp(CIA_Orton::CIA_H2_H2_e(It0, Inu0, nu_invcm, logT)) +
            (1 - mix) * exp(CIA_Orton::CIA_H2_H2_n(It0, Inu0, nu_invcm, logT)));
    h2he = amagat_h2 * amagat_he *
           (mix * exp(CIA_Orton::CIA_H2_He_e(It0, Inu0, nu_invcm, logT)) +
            (1 - mix) * exp(CIA_Orton::CIA_H2_He_n(It0, Inu0, nu_invcm, logT)));
    h2ch4 =
        amagat_h2 * amagat_ch4 *
        (mix * exp(CIA_Orton::CIA_H2_CH4_e(It0, Inu0, nu_invcm, logT)) +
         (1 - mix) * exp(CIA_Orton::CIA_H2_CH4_n(It0, Inu0, nu_invcm, logT)));
  }

  if (nu_scale != 1.0) {
    h2h2 *= nu_scale;
    h2he *= nu_scale;
    h2ch4 *= nu_scale;
  }

  return h2h2 + h2he + h2ch4;
}
