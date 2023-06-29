/** @file hydrogen.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Jul 10, 2021 14:55:03 PDT
 * @bug No known bugs.
 */

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include "molecules.hpp"

#define SQUARE(V) (V * V)
#define JMAX 20

static double const h_cgs = 6.62606957e-27;      // Planck constant
static double const c_cgs = 2.99792458e+10;      // speed of light
static double const kBoltz_cgs = 1.3806504e-16;  // erg/K
static double const FACT_cgs = h_cgs * c_cgs / kBoltz_cgs;
static double const R_cgs =
    8.314462e+07;  // universal gas constant in erg/K/mole

double Hydrogen::FJ_[JMAX];
double Hydrogen::DJ_[JMAX];
double Hydrogen::nist_shomate1_[7] = {33.066178, -11.363417, 11.432816,
                                      -2.772874, -0.158558,  -9.980797,
                                      172.707974};
double Hydrogen::nist_shomate2_[7] = {18.563083, 12.257357, -2.859786, 0.268238,
                                      1.977990,  -1.147438, 156.288133};
double Hydrogen::nist_shomate3_[7] = {43.413560, -4.293079,  1.272428,
                                      -0.096876, -20.533862, -38.515158,
                                      162.081354};

Hydrogen::Hydrogen() {
  double B0 = 59.322, D0 = 4.71e-02;

  for (int J = 0; J < JMAX; J++) {
    DJ_[J] = ((J / 2) * 2 == J) ? 2 * J + 1 : 3 * (2 * J + 1);
    FJ_[J] = B0 * J * (J + 1.) - D0 * SQUARE(J * (J + 1.));
  }
}

double Hydrogen::fpara_equil(double T) {
  // compute partition functions
  double Zpara = 0, Zortho = 0;

  for (size_t J = 0; J < JMAX; J += 2)
    Zpara += DJ_[J] * exp(-FJ_[J] * FACT_cgs / T);

  for (size_t J = 1; J < JMAX; J += 2)
    Zortho += DJ_[J] * exp(-FJ_[J] * FACT_cgs / T);

  return Zpara / (Zpara + Zortho);
}

void Hydrogen::get_cp_(double *cp_para, double *cp_equil, double *cp_norm,
                       double *cp_ortho, double T) {
  double FORTH = 0.75;
  double FPARA = 0.25;

  double ZEQ = 0.0;
  double SEQ1 = 0.0;
  double SEQ2 = 0.0;

  double ZORTHO = 0.0;
  double SO1 = 0.0;
  double SO2 = 0.0;

  double ZPARA = 0.0;
  double SP1 = 0.0;
  double SP2 = 0.0;

  for (int JQ = 9; JQ >= 0; JQ--) {
    double X = (FACT_cgs / T) * FJ_[JQ];
    double P = DJ_[JQ] * exp(-X);

    ZEQ += P;
    SEQ1 += DJ_[JQ] * SQUARE(X) * exp(-X);
    SEQ2 += DJ_[JQ] * X * exp(-X);

    if ((JQ / 2) * 2 == JQ) {
      ZPARA += P;
      SP1 += DJ_[JQ] * SQUARE(X) * exp(-X);
      SP2 += DJ_[JQ] * X * exp(-X);
    } else {
      ZORTHO += P;
      SO1 += DJ_[JQ] * SQUARE(X) * exp(-X);
      SO2 += DJ_[JQ] * X * exp(-X);
    }
  }

  *cp_equil = 5.0 / 2.0 + SEQ1 / ZEQ - SQUARE(SEQ2 / ZEQ);
  *cp_equil = (*cp_equil) * R_cgs * 1.E-7;

  *cp_ortho = 5.0 / 2.0 + SO1 / ZORTHO - SQUARE(SO2 / ZORTHO);
  *cp_ortho = (*cp_ortho) * R_cgs * 1.E-7;

  *cp_para = 5.0 / 2.0 + SP1 / ZPARA - SQUARE(SP2 / ZPARA);
  *cp_para = (*cp_para) * R_cgs * 1.E-7;

  *cp_norm = FORTH * (*cp_ortho) + FPARA * (*cp_para);
}

double Hydrogen::cp_para(double T) {
  double para;
  double equil;
  double norm;
  double ortho;
  get_cp_(&para, &equil, &norm, &ortho, T);
  return para;
}

double Hydrogen::cp_ortho(double T) {
  double para;
  double equil;
  double norm;
  double ortho;
  get_cp_(&para, &equil, &norm, &ortho, T);
  return ortho;
}

double Hydrogen::cp_norm(double T) {
  double para;
  double equil;
  double norm;
  double ortho;
  get_cp_(&para, &equil, &norm, &ortho, T);
  return norm;
}

double Hydrogen::cp_equil(double T) {
  double para;
  double equil;
  double norm;
  double ortho;
  get_cp_(&para, &equil, &norm, &ortho, T);
  return equil;
}

double Hydrogen::cp_nist(double T) {
  double *pdata;
  std::stringstream msg;
  T = std::min(std::max(298., T), 6000.);
  // if (T < 298. || T > 6000.) {
  //   msg << "ERROR: Temperature out of range in Hydrogen::cp_nist" <<
  //   std::endl; throw std::runtime_error(msg.str().c_str());
  // }

  if (T < 1000.)
    pdata = nist_shomate1_;
  else if (T >= 1000. && T < 2500.)
    pdata = nist_shomate2_;
  else
    pdata = nist_shomate3_;

  double result;
  T /= 1.E3;
  result = pdata[0] + pdata[1] * T + pdata[2] * T * T + pdata[3] * T * T * T +
           pdata[4] / (T * T);
  return result;
}

#undef SQUARE
#undef JMAX

Hydrogen aH2;
