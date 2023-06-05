#ifndef ABSORPTION_FUNCTIONS_HPP_
#define ABSORPTION_FUNCTIONS_HPP_

double absorption_coefficient_CIA(double freq, double P, double T, double XH2, double XHe,
       double XCH4, double mix);

double attenuation_NH3_Hanley(double freq, double P, double P_idl, double T,
       double XH2, double XHe, double XNH3, double XH2O = 0, double power = 0.);
double attenuation_NH3_Devaraj(double freq, double P, double P_idl, double T,
       double XH2, double XHe, double XNH3,double XH2O = 0, int version = 0);
double attenuation_NH3_Bellotti(double freq, double P, double P_idl, double T,
       double XH2, double xHe, double XNH3, double XH2O = 0);
double attenuation_NH3_Bellotti_switch(double freq, double P, double P_idl, double T,
       double XH2, double xHe, double XNH3, double XH2O = 0);
double attenuation_NH3_radtran(double freq, double P, double T, double XH2,
       double XHe, double XNH3);

double absorption_coefficient_PH3_radtran(double freq, double P, double T, double XH2, double XHe,
       double XPH3);

double absorption_coefficient_PH3_Hoffman(double freq, double P, double T, double XH2, double XHe,
       double XPH3);

double absorption_coefficient_H2S(double freq, double P, double T, double XH2, double XHe,
       double XH2S);

double attenuation_H2O_deBoer(double freq, double P, double T, double XH2, double XHe,
       double XH2O);
double attenuation_H2O_Waters(double freq, double P, double T, double XH2, double XHe,
       double XH2O);
double attenuation_H2O_Goodman(double freq, double P, double T, double XH2, double XHe,
       double XH2O);
double attenuation_H2O_Karpowicz(double freq, double P_idl, double T, double XH2,
       double XHe, double XH2O, double scale);

double absorption_coefficient_cloud(double freq, double P, double T, double rho_NH3_H2O,
  double rho_H2O, double rho_NH4SH, double rho_NH3, double cfliq, double cfwice, double cfaice);

double attenuation_freefree_Reference(double freq_GHz, double P_bar, double T);
double attenuation_freefree_Chengli(double freq_GHz, double P_bar, double T);
double attenuation_appleton_hartree_nomag(double freq_GHz, double P_bar, double T, double ne);

#endif
