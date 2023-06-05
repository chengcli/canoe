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
inline double SQUARE(double x) { return x*x; }
inline double CUBE(double x) { return x*x*x; }
inline double atm_from_bar(double P) { return 0.9869232667160128 * P; }
inline double bar_from_atm(double P) { return 1.013250 * P; }
static double const pi           = 3.14159265358979;
static double const kBoltz_mks   = 1.3806504e-23;      // J/K
static double invcm_from_dB_per_km = 2.3059e-6;    // dB/km to cm^-1
// end

double attenuation_H2O_Karpowicz(double freq, double P_idl, double T,
       double XH2, double XHe, double XH2O, double scale)
{
   double const C_H2Ob = 4.36510480961e-7; // (Mbar*GHz)^-2 (km)^-1
   double const C_H2Oa = 3.1e-7; // (Mbar*GHz)^-2 (km)^-1 (latest update: Nov 2012)
   double C_H2O = C_H2Oa*(1. - scale) + C_H2Ob*scale;

   double const Cprime_H2Ob = 2.10003048186e-26;
   double const Cprime_H2Oa = 0; // (latest update: Nov 2012)
   double Cprime_H2O = Cprime_H2Oa*(1. - scale) + Cprime_H2Ob*scale;

   double const C_H2 = 5.07722009423e-11; // (Mbar*GHz)^-2 (km)^-1
   double const C_He = 1.03562010226e-10; // (Mbar*GHz)^-2 (km)^-1
   double const x_ctmb = 13.3619799812;
   double const x_ctma = 12;  // (latest update: Nov 2012)
   double x_ctm = x_ctma*(1. - scale) + x_ctmb*scale;

   double const n_ctm = 6.76418487001;
   double const xprime_ctm = 0.0435525417274;

   // GHz
   static double const line[] = {
       22.2351,      183.3101,      321.2256,      325.1529,
      380.1974,      439.1508,      443.0183,      448.0011,
      470.8890,      474.6891,      488.4911,      556.9360,
      620.7008,      752.0332,      916.1712 };

   // Np-Hz/cm^2 (corigendum says this should be cm^2-Hz/molecule)
   static double const intensity[] = {
      0.1314e-13,      0.2279e-11,      0.8058e-13,      0.2701e-11,
      0.2444e-10,      0.2185e-11,      0.4637e-12,      0.2568e-10,
      0.8392e-12,      0.3272e-11,      0.6676e-12,      0.1535e-8,
      0.1711e-10,      0.1014e-8,       0.4238e-10 };

   // GHz/bar
   static double const linewidth_H2O[] = {
      13.49 ,      14.66 ,      10.57 ,      13.81 ,
      14.54 ,       9.715,       7.88 ,      12.75 ,
       9.83 ,      10.95 ,      13.13 ,      14.05 ,
      11.836,      12.53 ,      12.75 };

   // GHz/bar
   static double const linewidth_H2[]={
      2.395,      2.400,      2.395,      2.395,
      2.390,      2.395,      2.395,      2.395,
      2.395,      2.395,      2.395,      2.395,
      2.395,      2.395,      2.395 };

   // GHz/bar
   static double const linewidth_He[]={
      0.67,      0.71,      0.67,      0.67,
      0.63,      0.67,      0.67,      0.67,
      0.67,      0.67,      0.67,      0.67,
      0.67,      0.67,      0.67 };

   static double const exp_H2O[]={
      0.61,      0.85,      0.54,      0.74,
      0.89,      0.62,      0.50,      0.67,
      0.65,      0.64,      0.72,      1.0 ,
      0.68,      0.84,      0.78};

   static double const temp_coeff[]={
      2.144,      0.668,      6.179,      1.541,
      1.048,      3.595,      5.048,      1.405,
      3.597,      2.379,      2.852,      0.159,
      2.391,      0.396,      1.441};

   static double const exp_H2[] = {
      0.900,      0.950,      0.900,      0.900,
      0.850,      0.900,      0.900,      0.900,
      0.900,      0.900,      0.900,      0.900,
      0.900,      0.900,      0.900 };

   static double const exp_He[] = {
      0.515,      0.490,      0.515,      0.490,
      0.540,      0.515,      0.515,      0.515,
      0.515,      0.515,      0.515,      0.515,
      0.515,      0.515,      0.515 };

   double const PH2O_idl =  P_idl * XH2O;
   double const PH2_idl  =  P_idl * XH2;
   double const PHe_idl  =  P_idl * XHe;

   double const TR = 300 / T;
   double const nH2O = 0.1 * PH2O_idl / (kBoltz_mks*T); // molecules/cm^3

   // for now, assume only one isotope
   double abs_lines=0.0;
   for (auto i=0; i<sizeof(line) / sizeof(double); ++i) {
      // should have units of GHz
      double gamma = PH2O_idl*pow(TR,exp_H2O[i])*linewidth_H2O[i] +
	 PH2_idl*pow(TR,exp_H2[i])*linewidth_H2[i] +
	 PHe_idl*pow(TR,exp_He[i])*linewidth_He[i];

      // 1/GHz
      double sum_nu_line = fabs(freq + line[i]);
      double diff_nu_line = fabs(freq - line[i]);
      double cutoff_term = 1.0/(SQUARE(750.0) + SQUARE(gamma));

      double plus_term = ( sum_nu_line < 750
			   ? 1.0/(SQUARE(sum_nu_line) + SQUARE(gamma)) - cutoff_term
			   : 0.0);

      double minus_term = ( diff_nu_line < 750
			    ? 1.0/(SQUARE(diff_nu_line) + SQUARE(gamma)) - cutoff_term
			    : 0.0);

      // no cutoff
      //plus_term = 1.0/(SQUARE(sum_nu_line) + SQUARE(gamma));
      //minus_term = 1.0/(SQUARE(diff_nu_line) + SQUARE(gamma));

      double lineshape = SQUARE(freq/line[i]) * gamma/pi * (plus_term + minus_term);

      abs_lines += intensity[i] * exp((1-TR)*temp_coeff[i]) * lineshape;
   }
   abs_lines *= pow(TR,2.5) * nH2O;
   abs_lines *= 1e-9; // Hz/GHz --> unitless

   double mbar_of_bar = 1.0e3;
   double abs_ctm=0.0;
   abs_ctm += C_H2O * SQUARE(mbar_of_bar*PH2O_idl)*pow(TR,x_ctm);
   abs_ctm += Cprime_H2O * pow(mbar_of_bar*PH2O_idl,n_ctm)*pow(TR,xprime_ctm);
   abs_ctm += (C_H2*PH2_idl + C_He*PHe_idl) * PH2O_idl * CUBE(TR) * SQUARE(mbar_of_bar);
   abs_ctm *= SQUARE(freq);
   abs_ctm *= 1e-5; // inv(km) -> inv(cm)

   double abs_tot = abs_ctm + abs_lines; // 1/cm

   #ifdef ABSORPTION_KLUDGE
   constexpr double P_top = 150;
   constexpr double P_bot = 4000;
   constexpr double scale_bot = 1.0;
   if (P_idl>P_bot)
      abs_tot *= scale_bot;
   else if (P_idl>P_top)
      abs_tot *= ((P_idl-P_top)*scale_bot + (P_bot-P_idl)) / (P_bot-P_top);
   #endif

   return abs_tot;
}

double attenuation_H2O_Goodman(double freq, double P, double T, double XH2,
       double XHe, double XH2O)
{
   double const PH2O = atm_from_bar(P * XH2O);
   double const PH2  = atm_from_bar(P * XH2);
   double const PHe  = atm_from_bar(P * XHe);

   double TR = 273.0/T;
   double F2 = freq*freq;

   double DNU1 = 9.88e-2*pow(TR,0.6667)*(0.81*PH2 + 0.35*PHe);
   double DNU2 = DNU1*DNU1;
   double DNUP = freq/29.97 + 0.74;
   double DNUM = freq/29.97 - 0.74;
   double DTAUA = PH2O*F2*pow(TR,4.3333) *
      (9.07e-9*(DNU1/(DNUM*DNUM + DNU2) + DNU1/(DNUP*DNUP + DNU2)) + 1.45e-7*DNU1);

   return DTAUA;
}

double attenuation_H2O_Waters(double freq, double P, double T, double XH2,
       double XHe, double XH2O)
{
   double const P_atm = atm_from_bar(P);
   double const PH2O = atm_from_bar(P) * XH2O;

   double F2 = freq*freq;

   double TR = 273.0/T;
   double DTAUA = 0.5596 * F2 * P_atm * PH2O * pow(TR,3.126);
   // higher-order terms:
   double X1 = 1.0 + 3.897*PH2O/P_atm;
   double DNU1 = 2.96 * X1 * pow((300.0/T),0.626) * P_atm;
   double X2 = pow((494.4019-F2),2) + 4*F2*DNU1*DNU1;

   DTAUA *= X1* ( 2.77e-8 + 7.18*exp(-644.0/T)/(T*X2) );

   return DTAUA;
}


double attenuation_H2O_deBoer(double freq, double P, double T, double XH2,
       double XHe, double XH2O)
{
   static double const fi[10] = { 22.23515,183.31012,323.0,325.1538,
                                  380.1968,390.0,436.0,438.0,442.0,448.0 };

   static double const ei[10] = { 644.0,196.0,1850.0,454.0,306.0,2199.0,1507.0,
                                  1070.0,1507.0,412.0 };

   static double const ai[10] = { 1.0,41.9,334.4,115.7,651.8,127.0,191.4,
                                  697.6,590.2,973.1 };

   static double const c1[10][3] = { {2.395,0.67,10.67}, {2.40,0.71,11.64},
                                     {2.395,0.67,9.59}, {2.395,0.67,11.99},
                                     {2.39,0.63,12.42}, {2.395,0.67,9.16},
                                     {2.395,0.67,6.32}, {2.395,0.67,8.34},
                                     {2.395,0.67,6.52}, {2.395,0.67,11.57} };

   static double const c2[10][3] = { {0.9,0.515,0.626}, {0.95,0.49,0.649},
                                     {0.9,0.515,0.42}, {0.9,0.515,0.619},
                                     {0.85,0.54,0.63}, {0.9,0.515,0.33},
                                     {0.9,0.515,0.29}, {0.9,0.515,0.36},
                                     {0.9,0.515,0.332}, {0.9,0.515,0.51} };

   double const PH2O = atm_from_bar(P * XH2O);
   double const PH2  = atm_from_bar(P * XH2);
   double const PHe  = atm_from_bar(P * XHe);

   double TR = 300.0/T;
   double F2 = freq*freq;

   // convert P(atm) to P(bars) and store in convenient array:
   double px[3];
   px[0] = bar_from_atm(PH2);
   px[1] = bar_from_atm(PHe);
   px[2] = bar_from_atm(PH2O);

   double sum = 0.0;
   for (int i=0; i<10; i++) {
      double gami = 0.0;
      for (int j=0; j<3; j++) {
         gami += c1[i][j]*px[j]*pow(TR,c2[i][j]);
      }
      sum += ai[i]*exp(-ei[i]/T)*gami /
         ( pow((fi[i]*fi[i]-F2),2) + 4*F2*gami*gami );
   }

   double DTAUA = 1444.5 * px[2] * pow(TR,3.5) * F2 * sum;
   DTAUA += 0.00339 * px[2] * pow(TR,3.1) *
      (0.81*px[0]+0.35*px[1]) * F2;
   // convert dB/km to 1/cm:
   DTAUA *= invcm_from_dB_per_km;

   return DTAUA;
}
